//! Full non-linear multi-compartmental Hovorka (2004) endocrine and metabolic model.
//!
//! # Model Architecture
//! The Hovorka model consists of coupled physiological subsystems:
//! 1. Subcutaneous insulin absorption (two transit compartments $S_1, S_2$).
//! 2. Plasma insulin kinetics ($I(t)$).
//! 3. Three remote insulin action compartments:
//!    - $x_1$: Glucose transport into tissues.
//!    - $x_2$: Glucose disposal/utilization in tissues.
//!    - $x_3$: Inhibition of Endogenous Glucose Production ($EGP$).
//! 4. Two-compartment glucose subsystem:
//!    - $Q_1$: Accessible (vascular/interstitial) glucose mass.
//!    - $Q_2$: Non-accessible (intracellular/peripheral) glucose mass.
//! 5. Renal glucose clearance ($F_R$) when plasma glucose exceeds the tubular threshold.

use crate::error::{GlycoError, Result};
use crate::types::{HormoneInfusion, PatientProfile, INSULIN_U_TO_PMOL, MMOL_L_TO_MG_DL};

/// Hovorka model state vector length (8 primary state variables).
pub const HOVORKA_STATE_DIM: usize = 8;

/// State variable indices within the state vector.
pub mod state_idx {
    pub const S1: usize = 0; // Subcutaneous insulin compartment 1 (pmol)
    pub const S2: usize = 1; // Subcutaneous insulin compartment 2 (pmol)
    pub const I: usize = 2; // Plasma insulin concentration (pmol/L)
    pub const X1: usize = 3; // Remote insulin action on transport (1/min)
    pub const X2: usize = 4; // Remote insulin action on disposal (1/min)
    pub const X3: usize = 5; // Remote insulin action on EGP (unitless)
    pub const Q1: usize = 6; // Accessible glucose mass (mmol)
    pub const Q2: usize = 7; // Non-accessible glucose mass (mmol)
}

/// Biological parameters characterizing patient-specific Hovorka model kinetics.
#[derive(Debug, Clone)]
pub struct HovorkaParams {
    /// Subcutaneous insulin absorption time constant $\tau_s$ (min).
    pub tau_s: f64,
    /// Plasma insulin elimination rate constant $k_e$ (1/min).
    pub k_e: f64,
    /// Deactivation rate constant for insulin action on transport $k_{a1}$ (1/min).
    pub k_a1: f64,
    /// Activation rate constant for insulin action on transport $k_{b1}$ (1/min per pmol/L).
    pub k_b1: f64,
    /// Deactivation rate constant for insulin action on disposal $k_{a2}$ (1/min).
    pub k_a2: f64,
    /// Activation rate constant for insulin action on disposal $k_{b2}$ (1/min per pmol/L).
    pub k_b2: f64,
    /// Deactivation rate constant for insulin action on EGP $k_{a3}$ (1/min).
    pub k_a3: f64,
    /// Activation rate constant for insulin action on EGP $k_{b3}$ (1/min per pmol/L).
    pub k_b3: f64,
    /// Inter-compartmental glucose transfer rate constant $k_{12}$ (1/min).
    pub k_12: f64,
    /// Inter-compartmental glucose transfer rate constant $k_{21}$ (1/min).
    pub k_21: f64,
    /// Non-insulin-dependent glucose flux $F_{01}^c$ (mmol/min).
    pub f01_c: f64,
    /// Basal endogenous glucose production $EGP_0$ (mmol/min).
    pub egp_0: f64,
    /// Renal glucose clearance threshold (mmol/L).
    pub renal_threshold: f64,
    /// Renal clearance rate constant $c_R$ (1/min).
    pub renal_clearance_rate: f64,
    /// Glucose distribution volume $V_G$ (L).
    pub v_g: f64,
    /// Insulin distribution volume $V_I$ (L).
    pub v_i: f64,
}

impl HovorkaParams {
    /// Construct default parameters from a patient profile.
    pub fn from_patient(profile: &PatientProfile) -> Self {
        let v_g = profile.total_glucose_volume_l();
        let v_i = profile.total_insulin_volume_l();
        let tau_s = profile.tau_insulin_min;

        // Insulin sensitivity scalers derived from patient's ISF
        // Standard reference ISF is ~45 mg/dL/U
        let isf_ratio = (profile.isf_mg_dl_per_u / 45.0).clamp(0.4, 2.5);

        Self {
            tau_s,
            k_e: 0.138, // elimination half-life ~5 min
            k_a1: 0.006,
            k_b1: 5.12e-6 * isf_ratio,
            k_a2: 0.060,
            k_b2: 8.22e-6 * isf_ratio,
            k_a3: 0.030,
            k_b3: 2.60e-6 * isf_ratio,
            k_12: 0.066,
            k_21: 0.055,
            f01_c: profile.f01_c,
            egp_0: profile.egp_0,
            renal_threshold: profile.renal_threshold_mmol_l,
            renal_clearance_rate: 0.003,
            v_g,
            v_i,
        }
    }
}

/// Hovorka non-linear differential equation system.
#[derive(Debug, Clone)]
pub struct HovorkaSystem {
    pub params: HovorkaParams,
}

impl HovorkaSystem {
    /// Create new Hovorka system instance.
    pub fn new(params: HovorkaParams) -> Self {
        Self { params }
    }

    /// Evaluate system state derivatives $\frac{dy}{dt} = f(t, y)$.
    ///
    /// # Arguments
    /// * `t`: Current time (minutes).
    /// * `y`: State vector of length 8 ($S_1, S_2, I, x_1, x_2, x_3, Q_1, Q_2$).
    /// * `infusion`: Current hormone infusion commands.
    /// * `gut_absorption_rate_mmol_min`: Carbohydrate rate of appearance $U_G(t)$ from digestion (mmol/min).
    /// * `exercise_uptake_mmol_min`: Additional muscle glucose uptake from physical exertion.
    /// * `glucagon_egp_boost_mmol_min`: Additional hepatic glucose production stimulated by glucagon.
    /// * `dydt`: Destination derivative vector.
    pub fn evaluate_derivatives(
        &self,
        _t: f64,
        y: &[f64],
        infusion: &HormoneInfusion,
        gut_absorption_rate_mmol_min: f64,
        exercise_uptake_mmol_min: f64,
        glucagon_egp_boost_mmol_min: f64,
        dydt: &mut [f64],
    ) -> Result<()> {
        if y.len() < HOVORKA_STATE_DIM || dydt.len() < HOVORKA_STATE_DIM {
            return Err(GlycoError::InvalidParameter {
                parameter: "state_dim".to_string(),
                value: y.len() as f64,
                reason: format!("State dimension must be at least {}", HOVORKA_STATE_DIM),
            });
        }

        let p = &self.params;

        // Unpack states
        let s1 = y[state_idx::S1].max(0.0);
        let s2 = y[state_idx::S2].max(0.0);
        let insulin_plasma = y[state_idx::I].max(0.0);
        let x1 = y[state_idx::X1].max(0.0);
        let x2 = y[state_idx::X2].max(0.0);
        let x3 = y[state_idx::X3].max(0.0);
        let q1 = y[state_idx::Q1].max(0.0);
        let q2 = y[state_idx::Q2].max(0.0);

        // Subcutaneous insulin basal infusion rate (pmol/min)
        let basal_pmol_min = (infusion.insulin_basal_u_per_hr * INSULIN_U_TO_PMOL) / 60.0;
        let u_insulin = basal_pmol_min;

        // 1. Subcutaneous insulin absorption kinetics
        // dS1/dt = u(t) - S1 / tau_s
        dydt[state_idx::S1] = u_insulin - (s1 / p.tau_s);
        // dS2/dt = (S1 - S2) / tau_s
        dydt[state_idx::S2] = (s1 - s2) / p.tau_s;

        // 2. Plasma insulin concentration
        // dI/dt = (S2 / tau_s) / V_I - k_e * I
        let insulin_appearance = (s2 / p.tau_s) / p.v_i;
        dydt[state_idx::I] = insulin_appearance - p.k_e * insulin_plasma;

        // 3. Remote insulin action kinetics
        // dx1/dt = -k_a1 * x1 + k_b1 * I
        dydt[state_idx::X1] = -p.k_a1 * x1 + p.k_b1 * insulin_plasma;
        // dx2/dt = -k_a2 * x2 + k_b2 * I
        dydt[state_idx::X2] = -p.k_a2 * x2 + p.k_b2 * insulin_plasma;
        // dx3/dt = -k_a3 * x3 + k_b3 * I
        dydt[state_idx::X3] = -p.k_a3 * x3 + p.k_b3 * insulin_plasma;

        // 4. Glucose kinetics
        let plasma_conc_mmol_l = q1 / p.v_g;

        // Renal clearance F_R
        let f_r = if plasma_conc_mmol_l > p.renal_threshold {
            p.renal_clearance_rate * (plasma_conc_mmol_l - p.renal_threshold) * p.v_g
        } else {
            0.0
        };

        // Non-insulin-dependent CNS glucose uptake F_01_c
        let f01 = if plasma_conc_mmol_l >= 4.5 {
            p.f01_c
        } else {
            p.f01_c * (plasma_conc_mmol_l / 4.5)
        };

        // Endogenous glucose production: EGP_0 * (1 - x3) + glucagon stimulation
        let egp = (p.egp_0 * (1.0 - x3).max(0.0) + glucagon_egp_boost_mmol_min).max(0.0);

        // Rate of appearance from gut
        let u_g = gut_absorption_rate_mmol_min;

        // Inter-compartmental glucose exchange flux (Q1 -> Q2): x1 * Q1 - k12 * Q2
        let inter_compartment_flux = (x1 * q1) - (p.k_12 * q2);

        // Total peripheral tissue disposal from non-accessible compartment Q2
        let disposal_rate = (x2 * q2) + exercise_uptake_mmol_min;

        // dQ1/dt = - F01 - inter_compartment_flux - F_R + U_G + EGP
        dydt[state_idx::Q1] = -f01 - inter_compartment_flux - f_r + u_g + egp;

        // dQ2/dt = inter_compartment_flux - disposal
        dydt[state_idx::Q2] = inter_compartment_flux - disposal_rate;

        Ok(())
    }

    /// Extract plasma glucose concentration in mg/dL from the state vector.
    pub fn plasma_glucose_mg_dl(&self, y: &[f64]) -> f64 {
        let q1 = y[state_idx::Q1].max(0.0);
        let mmol_l = q1 / self.params.v_g;
        mmol_l * MMOL_L_TO_MG_DL
    }

    /// Extract plasma insulin concentration in mU/L from the state vector (1 mU/L = 6.0 pmol/L).
    pub fn plasma_insulin_mu_l(&self, y: &[f64]) -> f64 {
        let pmol_l = y[state_idx::I].max(0.0);
        pmol_l / 6.0
    }
}
