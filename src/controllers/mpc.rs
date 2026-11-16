//! Constrained Model Predictive Control (MPC) with asymmetric cost functions for artificial pancreas delivery.
//!
//! # Objective Function
//! $$ J(u) = \sum_{k=1}^{N_p} \left[ q_k \cdot (G_k - G_{target})^2 + \Psi(G_k) \right] + \sum_{j=0}^{N_c-1} r_j \cdot (\Delta u_j)^2 $$
//! where $\Psi(G_k)$ is a steep asymmetric penalty for hypoglycemia:
//! $$ \Psi(G_k) = \begin{cases} \alpha_{hypo} \cdot (70 - G_k)^2 & \text{if } G_k < 70 \text{ mg/dL} \\ 0 & \text{otherwise} \end{cases} $$

use crate::types::{HormoneInfusion, PatientProfile};

/// MPC tuning parameters.
#[derive(Debug, Clone)]
pub struct MpcConfig {
    /// Prediction horizon length $N_p$ (steps, e.g. 24 steps of 5 min = 120 min).
    pub prediction_horizon: usize,
    /// Control horizon length $N_c$ (steps, e.g. 6 steps = 30 min).
    pub control_horizon: usize,
    /// Sample step size $\Delta t$ (minutes, e.g. 5.0 min).
    pub dt_min: f64,
    /// Hyperglycemia tracking error weight $Q$.
    pub q_tracking: f64,
    /// Asymmetric hypoglycemia penalty weight $\alpha_{hypo}$.
    pub alpha_hypo: f64,
    /// Slew rate penalty weight $R$ penalizing rapid changes in infusion rate $\Delta u$.
    pub r_slew: f64,
    /// Maximum allowed basal delivery rate multiplier relative to patient basal (e.g. 3.5x).
    pub max_basal_multiplier: f64,
    /// Hypoglycemia threshold setpoint (mg/dL).
    pub hypo_threshold_mg_dl: f64,
}

impl Default for MpcConfig {
    fn default() -> Self {
        Self {
            prediction_horizon: 24, // 120 min
            control_horizon: 6,     // 30 min
            dt_min: 5.0,
            q_tracking: 1.0,
            alpha_hypo: 18.0, // steep 18x penalty for hypoglycemia
            r_slew: 2.5,
            max_basal_multiplier: 3.5,
            hypo_threshold_mg_dl: 70.0,
        }
    }
}

/// Constrained Model Predictive Controller for artificial pancreas automated delivery.
#[derive(Debug, Clone)]
pub struct ModelPredictiveController {
    pub config: MpcConfig,
    pub patient: PatientProfile,
    last_commanded_u: f64,
}

impl ModelPredictiveController {
    /// Create new MPC controller for a given patient.
    pub fn new(config: MpcConfig, patient: PatientProfile) -> Self {
        let basal = patient.basal_rate_u_per_hr;
        Self {
            config,
            patient,
            last_commanded_u: basal,
        }
    }

    /// Compute optimal insulin and glucagon infusion commands given current CGM reading and trend.
    ///
    /// # Arguments
    /// * `current_glucose_mg_dl`: Latest CGM reading.
    /// * `glucose_trend_mg_dl_min`: Estimated rate of change $dG/dt$ (mg/dL/min).
    /// * `insulin_on_board_u`: Current active insulin on board ($IOB$).
    pub fn compute_action(
        &mut self,
        current_glucose_mg_dl: f64,
        glucose_trend_mg_dl_min: f64,
        insulin_on_board_u: f64,
    ) -> HormoneInfusion {
        let cfg = &self.config;
        let basal = self.patient.basal_rate_u_per_hr;
        let target = self.patient.target_glucose_mg_dl;
        let isf = self.patient.isf_mg_dl_per_u;

        let u_min = 0.0;
        let u_max = basal * cfg.max_basal_multiplier;

        // Optimization over control candidate inputs u in [u_min, u_max]
        // Discrete active-set grid search over candidate increments
        let candidates = 41;
        let mut best_u = basal;
        let mut min_cost = f64::MAX;

        let step_size = (u_max - u_min) / (candidates - 1) as f64;

        for c in 0..candidates {
            let u_cand = u_min + c as f64 * step_size;
            let delta_u = u_cand - self.last_commanded_u;

            // Slew penalty: r * (delta_u)^2
            let mut cost = cfg.r_slew * delta_u * delta_u;

            // Forward simulate linear impulse-response prediction model over N_p horizon
            for k in 1..=cfg.prediction_horizon {
                let t_k = k as f64 * cfg.dt_min;

                // Estimated glucose trajectory:
                // G_k = G_0 + t_k * (dG/dt) - ISF * (additional_insulin_effect)
                let delta_insulin_rate = u_cand - basal; // U/hr
                let delivered_extra_u = (delta_insulin_rate / 60.0) * t_k;

                // Insulin pharmacodynamic fraction decayed by DIA
                let dia = self.patient.dia_minutes;
                let active_fraction = (t_k / dia).clamp(0.0, 1.0);
                let pd_effect = delivered_extra_u * active_fraction * isf;

                // Damping trend projection as t increases
                let trend_decay = (-t_k / 45.0).exp();
                let projected_g = current_glucose_mg_dl
                    + (glucose_trend_mg_dl_min * 45.0 * (1.0 - trend_decay))
                    - pd_effect
                    - (insulin_on_board_u * isf * (1.0 - (-t_k / dia).exp()) * 0.4);

                // Tracking cost
                let err = projected_g - target;
                cost += cfg.q_tracking * err * err;

                // Asymmetric hypoglycemia barrier penalty
                if projected_g < cfg.hypo_threshold_mg_dl {
                    let hypo_deficit = cfg.hypo_threshold_mg_dl - projected_g;
                    cost += cfg.alpha_hypo * hypo_deficit * hypo_deficit;
                }
            }

            if cost < min_cost {
                min_cost = cost;
                best_u = u_cand;
            }
        }

        self.last_commanded_u = best_u;

        // Glucagon micro-dose logic: trigger rescue bolus if predicted nadir is critically low (< 60 mg/dL)
        let projected_nadir = current_glucose_mg_dl + glucose_trend_mg_dl_min * 20.0;
        let glucagon_bolus = if projected_nadir < 60.0 && glucose_trend_mg_dl_min < -0.5 {
            let deficit = 70.0 - projected_nadir;
            // 1 microgram of glucagon raises blood glucose by ~1.5 to 2.5 mg/dL in typical adults
            (deficit * 1.5).clamp(5.0, 35.0)
        } else {
            0.0
        };

        HormoneInfusion {
            insulin_basal_u_per_hr: best_u,
            insulin_bolus_u: 0.0,
            glucagon_basal_ug_per_hr: 0.0,
            glucagon_bolus_ug: glucagon_bolus,
        }
    }
}
