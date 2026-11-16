//! Glucagon pharmacokinetics, hepatic glycogenolysis saturation, receptor desensitization, and neuroendocrine counter-regulation.

use crate::types::HormoneInfusion;

/// Glucagon physiological parameters.
#[derive(Debug, Clone)]
pub struct GlucagonParams {
    /// Subcutaneous glucagon absorption time constant (min).
    pub tau_glucagon_min: f64,
    /// Glucagon elimination half-life rate constant (1/min, $t_{1/2} \approx 6$ min).
    pub k_clearance: f64,
    /// Basal plasma glucagon concentration (pg/mL, typically 50-80 pg/mL).
    pub basal_plasma_conc: f64,
    /// Apparent distribution volume for glucagon (L/kg).
    pub distribution_volume_l_per_kg: f64,
    /// Maximum hepatic glucose output rate stimulated by glucagon ($E_{max}$, mmol/min).
    pub max_hepatic_output_mmol_min: f64,
    /// Concentration of glucagon producing 50% of maximum effect ($EC_{50}$, pg/mL).
    pub ec50_pg_ml: f64,
    /// Hill coefficient for hepatic glycogenolysis activation.
    pub hill_coefficient: f64,
    /// Glycogen pool maximum capacity (mmol of glucose equivalents, ~500 mmol in liver).
    pub max_glycogen_pool_mmol: f64,
    /// Rate of glycogen replenishment from circulating glucose (mmol/min).
    pub glycogen_synthesis_rate: f64,
    /// Receptor desensitization / tachyphylaxis time constant (min).
    pub tau_desensitization_min: f64,
}

impl Default for GlucagonParams {
    fn default() -> Self {
        Self {
            tau_glucagon_min: 18.0,
            k_clearance: 0.12, // t_1/2 ~ 5.8 min
            basal_plasma_conc: 65.0,
            distribution_volume_l_per_kg: 0.18,
            max_hepatic_output_mmol_min: 1.85,
            ec50_pg_ml: 120.0,
            hill_coefficient: 2.2,
            max_glycogen_pool_mmol: 500.0,
            glycogen_synthesis_rate: 0.15,
            tau_desensitization_min: 180.0,
        }
    }
}

/// Dynamic state of the glucagon counter-regulatory subsystem.
#[derive(Debug, Clone)]
pub struct GlucagonState {
    /// Subcutaneous glucagon depot 1 ($\mu$g).
    pub s1_subcut_ug: f64,
    /// Subcutaneous glucagon depot 2 ($\mu$g).
    pub s2_subcut_ug: f64,
    /// Plasma glucagon concentration (pg/mL).
    pub plasma_glucagon_pg_ml: f64,
    /// Current available liver glycogen pool (fraction of capacity 0.0 to 1.0).
    pub glycogen_fraction: f64,
    /// Receptor sensitivity multiplier (1.0 = fully sensitized, drops during prolonged high dose).
    pub receptor_sensitivity: f64,
    /// Cumulative epinephrine/stress hormone surge multiplier (1.0 = resting, > 1.0 = acute stress).
    pub stress_multiplier: f64,
}

impl GlucagonState {
    /// Initial baseline state.
    pub fn new(params: &GlucagonParams) -> Self {
        Self {
            s1_subcut_ug: 0.0,
            s2_subcut_ug: 0.0,
            plasma_glucagon_pg_ml: params.basal_plasma_conc,
            glycogen_fraction: 0.85,
            receptor_sensitivity: 1.0,
            stress_multiplier: 1.0,
        }
    }
}

/// Glucagon and neuroendocrine counter-regulation controller and simulation module.
#[derive(Debug, Clone)]
pub struct CounterRegulationSystem {
    pub params: GlucagonParams,
    pub state: GlucagonState,
}

impl CounterRegulationSystem {
    /// Construct a new counter-regulation system.
    pub fn new(params: GlucagonParams) -> Self {
        let state = GlucagonState::new(&params);
        Self { params, state }
    }

    /// Advance counter-regulatory state across time interval `dt` (min).
    ///
    /// Returns the net hepatic glucose production boost (mmol/min) to be added to EGP.
    pub fn step(
        &mut self,
        dt: f64,
        infusion: &HormoneInfusion,
        body_weight_kg: f64,
        plasma_glucose_mg_dl: f64,
        plasma_insulin_mu_l: f64,
    ) -> f64 {
        let p = &self.params;
        let s = &mut self.state;

        // Subcutaneous delivery input (ug/min)
        let basal_ug_min = infusion.glucagon_basal_ug_per_hr / 60.0;
        // Deliver discrete bolus directly into subcutaneous depot
        s.s1_subcut_ug += infusion.glucagon_bolus_ug;
        let u_glucagon = basal_ug_min;

        // Two-compartment subcutaneous absorption
        let ds1 = (u_glucagon - (s.s1_subcut_ug / p.tau_glucagon_min)) * dt;
        s.s1_subcut_ug = (s.s1_subcut_ug + ds1).max(0.0);

        let ds2 = ((s.s1_subcut_ug - s.s2_subcut_ug) / p.tau_glucagon_min) * dt;
        s.s2_subcut_ug = (s.s2_subcut_ug + ds2).max(0.0);

        // Systemic appearance (1 ug = 1,000,000 pg)
        let volume_l = p.distribution_volume_l_per_kg * body_weight_kg;
        let volume_ml = volume_l * 1000.0;
        let rate_of_appearance_pg_ml_min = (s.s2_subcut_ug / p.tau_glucagon_min * 1e6) / volume_ml;

        // Endogenous alpha-cell secretion in response to hypoglycemia
        // Stimulated when plasma glucose drops below 75 mg/dL
        let endogenous_alpha_secretion = if plasma_glucose_mg_dl < 75.0 {
            let deficit = 75.0 - plasma_glucose_mg_dl;
            // Suppressed by high ambient insulin
            let insulin_suppression = 1.0 / (1.0 + (plasma_insulin_mu_l / 15.0).powi(2));
            (deficit * 2.5) * insulin_suppression
        } else {
            0.0
        };

        // Update plasma glucagon
        let d_plasma = (rate_of_appearance_pg_ml_min + endogenous_alpha_secretion
            - p.k_clearance * (s.plasma_glucagon_pg_ml - p.basal_plasma_conc))
            * dt;
        s.plasma_glucagon_pg_ml =
            (s.plasma_glucagon_pg_ml + d_plasma).max(p.basal_plasma_conc * 0.5);

        // Receptor desensitization under elevated glucagon
        if s.plasma_glucagon_pg_ml > p.basal_plasma_conc * 2.0 {
            let d_sens = -(s.plasma_glucagon_pg_ml / (p.ec50_pg_ml * 5.0))
                * (dt / p.tau_desensitization_min);
            s.receptor_sensitivity = (s.receptor_sensitivity + d_sens).clamp(0.20, 1.0);
        } else {
            // Gradual resensitization
            let d_sens = 0.05 * (1.0 - s.receptor_sensitivity) * (dt / p.tau_desensitization_min);
            s.receptor_sensitivity = (s.receptor_sensitivity + d_sens).clamp(0.20, 1.0);
        }

        // Hepatic glucose production via Hill equation
        let delta_conc = (s.plasma_glucagon_pg_ml - p.basal_plasma_conc).max(0.0);
        let hill_term = if delta_conc > 0.0 {
            delta_conc.powf(p.hill_coefficient)
                / (p.ec50_pg_ml.powf(p.hill_coefficient) + delta_conc.powf(p.hill_coefficient))
        } else {
            0.0
        };

        // Glycogen availability gate (depleted glycogen limits release)
        let glycogen_factor = s.glycogen_fraction.clamp(0.0, 1.0);

        // Insulin suppression of glucagon action: high insulin blunts hepatic responsiveness
        let insulin_brake = 1.0 / (1.0 + (plasma_insulin_mu_l / 25.0));

        // Net hepatic glucose output boost (mmol/min)
        let hepatic_boost = p.max_hepatic_output_mmol_min
            * hill_term
            * s.receptor_sensitivity
            * glycogen_factor
            * insulin_brake
            * s.stress_multiplier;

        // Deplete liver glycogen by hepatic glucose release, replenished continuously
        let glycogen_depletion = (hepatic_boost * dt) / p.max_glycogen_pool_mmol;
        let glycogen_replenishment = (p.glycogen_synthesis_rate * dt) / p.max_glycogen_pool_mmol;
        s.glycogen_fraction =
            (s.glycogen_fraction - glycogen_depletion + glycogen_replenishment).clamp(0.05, 1.0);

        hepatic_boost
    }
}
