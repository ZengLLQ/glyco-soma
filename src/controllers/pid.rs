//! Dual-hormone Proportional-Integral-Derivative (PID) controller with anti-windup clamping and noise filtering.

use crate::types::{HormoneInfusion, PatientProfile};

/// PID controller configuration parameters.
#[derive(Debug, Clone)]
pub struct PidConfig {
    /// Proportional gain $K_p$ (U/hr per mg/dL error).
    pub kp: f64,
    /// Integral gain $K_i$ (U/hr per mg/dL * min).
    pub ki: f64,
    /// Derivative gain $K_d$ (U/hr per mg/dL / min).
    pub kd: f64,
    /// Derivative low-pass filter time constant $\tau_d$ (min, typically 2.0 to 5.0 min).
    pub tau_filter_min: f64,
    /// Anti-windup maximum integral accumulator magnitude.
    pub max_integral_error: f64,
    /// Maximum allowed total basal delivery rate (U/hr).
    pub max_delivery_u_per_hr: f64,
    /// Minimum allowed total basal delivery rate (U/hr, typically 0.0).
    pub min_delivery_u_per_hr: f64,
    /// Glucagon proportional gain $K_{p,gluc}$ ($\mu$g/hr per mg/dL below target).
    pub kp_glucagon: f64,
}

impl Default for PidConfig {
    fn default() -> Self {
        Self {
            kp: 0.015,
            ki: 0.00015,
            kd: 0.12,
            tau_filter_min: 3.5,
            max_integral_error: 450.0,
            max_delivery_u_per_hr: 4.5,
            min_delivery_u_per_hr: 0.0,
            kp_glucagon: 0.45,
        }
    }
}

/// Dual-hormone PID feedback controller.
#[derive(Debug, Clone)]
pub struct PidController {
    pub config: PidConfig,
    pub patient: PatientProfile,
    integral_error: f64,
    last_error: f64,
    filtered_derivative: f64,
    last_time_min: f64,
}

impl PidController {
    /// Construct a new PID controller for a given patient.
    pub fn new(config: PidConfig, patient: PatientProfile) -> Self {
        Self {
            config,
            patient,
            integral_error: 0.0,
            last_error: 0.0,
            filtered_derivative: 0.0,
            last_time_min: 0.0,
        }
    }

    /// Reset internal integral and derivative states.
    pub fn reset(&mut self) {
        self.integral_error = 0.0;
        self.last_error = 0.0;
        self.filtered_derivative = 0.0;
        self.last_time_min = 0.0;
    }

    /// Compute PID action given current glucose measurement at time `t`.
    pub fn step(&mut self, t: f64, glucose_mg_dl: f64) -> HormoneInfusion {
        let cfg = &self.config;
        let dt = (t - self.last_time_min).clamp(0.5, 30.0);
        self.last_time_min = t;

        let target = self.patient.target_glucose_mg_dl;
        let error = glucose_mg_dl - target;

        // 1. Proportional term
        let p_term = cfg.kp * error;

        // 2. Integral term with anti-windup clamping
        // Suspend integration if in hypoglycemia or saturated
        if glucose_mg_dl > 70.0 {
            self.integral_error += error * dt;
            self.integral_error = self
                .integral_error
                .clamp(-cfg.max_integral_error, cfg.max_integral_error);
        } else {
            // Rapid decay of integral in hypoglycemia to prevent stacking
            self.integral_error *= 0.5;
        }
        let i_term = cfg.ki * self.integral_error;

        // 3. Derivative term with first-order low-pass noise filter
        // dE_filt / dt = (dE/dt - E_filt) / tau_filter
        let raw_derivative = (error - self.last_error) / dt;
        self.last_error = error;

        let filter_alpha = dt / (dt + cfg.tau_filter_min);
        self.filtered_derivative =
            (1.0 - filter_alpha) * self.filtered_derivative + filter_alpha * raw_derivative;
        let d_term = cfg.kd * self.filtered_derivative;

        // Net control delta relative to patient baseline basal
        let nominal_basal = self.patient.basal_rate_u_per_hr;
        let commanded_u = nominal_basal + p_term + i_term + d_term;
        let final_insulin = commanded_u.clamp(cfg.min_delivery_u_per_hr, cfg.max_delivery_u_per_hr);

        // Dual-hormone glucagon action:
        // Activate glucagon when glucose is below 75 mg/dL with negative trend
        let glucagon_bolus = if glucose_mg_dl < 75.0 && self.filtered_derivative < 0.0 {
            let hypo_deficit = 80.0 - glucose_mg_dl;
            let bolus_ug = hypo_deficit * cfg.kp_glucagon * 2.0;
            bolus_ug.clamp(0.0, 30.0)
        } else {
            0.0
        };

        HormoneInfusion {
            insulin_basal_u_per_hr: final_insulin,
            insulin_bolus_u: 0.0,
            glucagon_basal_ug_per_hr: 0.0,
            glucagon_bolus_ug: glucagon_bolus,
        }
    }
}
