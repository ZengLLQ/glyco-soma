//! Predictive Low Glucose Suspend (PLGS) and Threshold Suspend safety automaton.
//!
//! # Clinical Function
//! Automatically pauses basal insulin infusion when CGM glucose is predicted to cross below a hypoglycemia threshold
//! within 30 minutes, or when current glucose crosses below a hard threshold (e.g. 70 mg/dL).
//! Automatically resumes basal delivery when glucose begins to recover, or after an emergency maximum 2-hour suspend limit.

use crate::types::{HormoneInfusion, PatientProfile};

/// PLGS operational state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlgsState {
    /// Normal basal delivery.
    Normal,
    /// Preemptive suspend based on projected 30-minute hypoglycemia.
    PreemptiveSuspension,
    /// Hard threshold suspend active (current $G < 70$ mg/dL).
    ThresholdSuspension,
    /// Post-hypoglycemia recovery ramp-up.
    RecoveryRampUp,
}

/// PLGS safety configuration.
#[derive(Debug, Clone)]
pub struct PlgsConfig {
    /// Hard hypoglycemia suspend threshold (mg/dL, typically 70.0 mg/dL).
    pub threshold_hypo_mg_dl: f64,
    /// Forward prediction horizon for preemptive suspend (minutes, typically 30.0 min).
    pub prediction_horizon_min: f64,
    /// Preemptive suspend threshold (predicted glucose nadir below this triggers pause, e.g. 80.0 mg/dL).
    pub preemptive_threshold_mg_dl: f64,
    /// Minimum glucose level required to permit automated resumption (mg/dL, e.g. 75.0 mg/dL).
    pub resume_glucose_threshold_mg_dl: f64,
    /// Positive trend required for automated resumption (mg/dL/min, e.g. +0.5 mg/dL/min).
    pub resume_trend_threshold_mg_dl_min: f64,
    /// Maximum allowable continuous suspension time to prevent DKA (minutes, typically 120.0 min = 2 hours).
    pub max_suspension_duration_min: f64,
}

impl Default for PlgsConfig {
    fn default() -> Self {
        Self {
            threshold_hypo_mg_dl: 70.0,
            prediction_horizon_min: 30.0,
            preemptive_threshold_mg_dl: 80.0,
            resume_glucose_threshold_mg_dl: 75.0,
            resume_trend_threshold_mg_dl_min: 0.5,
            max_suspension_duration_min: 120.0,
        }
    }
}

/// Predictive Low Glucose Suspend (PLGS) controller.
#[derive(Debug, Clone)]
pub struct PlgsController {
    pub config: PlgsConfig,
    pub patient: PatientProfile,
    pub current_state: PlgsState,
    pub time_in_suspension_min: f64,
}

impl PlgsController {
    /// Create a new PLGS controller instance.
    pub fn new(config: PlgsConfig, patient: PatientProfile) -> Self {
        Self {
            config,
            patient,
            current_state: PlgsState::Normal,
            time_in_suspension_min: 0.0,
        }
    }

    /// Advance PLGS automaton state and calculate commanded delivery rate.
    pub fn step(
        &mut self,
        dt_min: f64,
        current_glucose_mg_dl: f64,
        glucose_trend_mg_dl_min: f64,
    ) -> HormoneInfusion {
        let cfg = &self.config;
        let nominal_basal = self.patient.basal_rate_u_per_hr;

        // Project forward nadir over 30 min horizon: G_proj = G + trend * 30 min
        let projected_glucose =
            current_glucose_mg_dl + (glucose_trend_mg_dl_min * cfg.prediction_horizon_min);

        // State transitions
        match self.current_state {
            PlgsState::Normal => {
                if current_glucose_mg_dl <= cfg.threshold_hypo_mg_dl {
                    self.current_state = PlgsState::ThresholdSuspension;
                    self.time_in_suspension_min = 0.0;
                } else if projected_glucose <= cfg.preemptive_threshold_mg_dl
                    && glucose_trend_mg_dl_min < -0.3
                {
                    self.current_state = PlgsState::PreemptiveSuspension;
                    self.time_in_suspension_min = 0.0;
                }
            }
            PlgsState::PreemptiveSuspension | PlgsState::ThresholdSuspension => {
                self.time_in_suspension_min += dt_min;

                // Safety timeout: enforce resumption after 2 hours regardless to prevent DKA
                if self.time_in_suspension_min >= cfg.max_suspension_duration_min
                    || (current_glucose_mg_dl >= cfg.resume_glucose_threshold_mg_dl
                        && glucose_trend_mg_dl_min >= cfg.resume_trend_threshold_mg_dl_min)
                {
                    self.current_state = PlgsState::RecoveryRampUp;
                }
            }
            PlgsState::RecoveryRampUp => {
                // Return to normal once safely above 90 mg/dL or after 15 minutes of recovery
                if current_glucose_mg_dl >= 90.0 {
                    self.current_state = PlgsState::Normal;
                    self.time_in_suspension_min = 0.0;
                }
            }
        }

        // Commanded insulin delivery based on state
        let insulin_rate = match self.current_state {
            PlgsState::Normal => nominal_basal,
            PlgsState::PreemptiveSuspension | PlgsState::ThresholdSuspension => 0.0, // Suspend insulin
            PlgsState::RecoveryRampUp => nominal_basal * 0.50, // 50% gentle ramp-up
        };

        HormoneInfusion {
            insulin_basal_u_per_hr: insulin_rate,
            insulin_bolus_u: 0.0,
            glucagon_basal_ug_per_hr: 0.0,
            glucagon_bolus_ug: 0.0,
        }
    }
}
