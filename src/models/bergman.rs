//! Classic Bergman Minimal Model for glucose-insulin interaction benchmarking.
//!
//! # Equations
//! $$ \frac{dG}{dt} = -(p_1 + X) G(t) + p_1 G_b + \frac{R_a(t)}{V_G} $$
//! $$ \frac{dX}{dt} = -p_2 X(t) + p_3 (I(t) - I_b) $$
//! $$ \frac{dI}{dt} = -n (I(t) - I_b) + \frac{u(t)}{V_I} $$

use crate::error::{GlycoError, Result};
use crate::types::PatientProfile;

/// Bergman minimal model parameters.
#[derive(Debug, Clone)]
pub struct BergmanParams {
    /// Glucose effectiveness $p_1$ (1/min).
    pub p1: f64,
    /// Insulin deactivation rate constant $p_2$ (1/min).
    pub p2: f64,
    /// Insulin sensitivity parameter $p_3$ (L / (mU * min^2)).
    pub p3: f64,
    /// Basal glucose concentration $G_b$ (mg/dL).
    pub g_b: f64,
    /// Basal insulin concentration $I_b$ (mU/L).
    pub i_b: f64,
    /// Insulin clearance rate constant $n$ (1/min).
    pub n: f64,
    /// Glucose distribution volume $V_G$ (dL).
    pub v_g: f64,
    /// Insulin distribution volume $V_I$ (L).
    pub v_i: f64,
}

impl BergmanParams {
    /// Create parameters matching patient profile.
    pub fn from_patient(profile: &PatientProfile) -> Self {
        let v_g_dl = profile.total_glucose_volume_l() * 10.0;
        let v_i_l = profile.total_insulin_volume_l();

        Self {
            p1: 0.028,
            p2: 0.025,
            p3: 0.000013 * (profile.isf_mg_dl_per_u / 45.0),
            g_b: profile.target_glucose_mg_dl,
            i_b: 11.0,
            n: 0.23,
            v_g: v_g_dl,
            v_i: v_i_l,
        }
    }
}

/// Bergman Minimal Model system.
#[derive(Debug, Clone)]
pub struct BergmanSystem {
    pub params: BergmanParams,
}

impl BergmanSystem {
    /// Create new Bergman system.
    pub fn new(params: BergmanParams) -> Self {
        Self { params }
    }

    /// Evaluate derivatives for state $[G, X, I]$:
    /// - $y[0]$: Glucose concentration $G$ (mg/dL)
    /// - $y[1]$: Remote insulin action $X$ (1/min)
    /// - $y[2]$: Plasma insulin $I$ (mU/L)
    pub fn evaluate_derivatives(
        &self,
        y: &[f64],
        insulin_infusion_u_per_hr: f64,
        glucose_appearance_mg_dl_min: f64,
        dydt: &mut [f64],
    ) -> Result<()> {
        if y.len() < 3 || dydt.len() < 3 {
            return Err(GlycoError::InvalidParameter {
                parameter: "y.len".to_string(),
                value: y.len() as f64,
                reason: "Bergman state requires at least 3 elements".to_string(),
            });
        }

        let p = &self.params;
        let g = y[0].max(10.0);
        let x = y[1];
        let i = y[2].max(0.0);

        // dG/dt = -(p1 + X)*G + p1*Gb + Ra
        dydt[0] = -(p.p1 + x) * g + p.p1 * p.g_b + glucose_appearance_mg_dl_min;

        // dX/dt = -p2*X + p3*(I - Ib)
        dydt[1] = -p.p2 * x + p.p3 * (i - p.i_b);

        // dI/dt = -n*(I - Ib) + u / VI
        let u_mu_min = (insulin_infusion_u_per_hr * 1000.0) / 60.0;
        dydt[2] = -p.n * (i - p.i_b) + (u_mu_min / p.v_i);

        Ok(())
    }

    /// Compute steady-state insulin infusion rate for basal glucose $G_b$.
    pub fn basal_infusion_rate_u_hr(&self) -> f64 {
        let u_mu_min = self.params.n * self.params.i_b * self.params.v_i;
        (u_mu_min * 60.0) / 1000.0
    }
}
