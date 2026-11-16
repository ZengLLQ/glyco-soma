//! Physical exertion metabolic kinetics, non-insulin mediated glucose uptake (NIMGU), and prolonged nocturnal insulin sensitivity.

use crate::types::ExerciseEvent;

/// Metabolic parameters governing physical exercise dynamics.
#[derive(Debug, Clone)]
pub struct ExerciseParams {
    /// Non-insulin-mediated glucose uptake multiplier per MET above resting (mmol/min/MET).
    pub nimgu_rate_per_met: f64,
    /// Muscle glycogen depletion constant per minute of high-intensity exertion.
    pub glycogen_depletion_rate: f64,
    /// Exercise-induced insulin sensitivity amplification factor ($S_{I,exercise} \approx 1.25$ to 2.0).
    pub sensitivity_boost_max: f64,
    /// Recovery decay time constant for acute post-exercise insulin sensitivity (minutes, $\tau \approx 120$ min).
    pub tau_acute_recovery_min: f64,
    /// Delayed nocturnal insulin sensitivity surge onset time after exercise (minutes, ~480 min = 8 hours).
    pub nocturnal_surge_onset_min: f64,
    /// Delayed nocturnal insulin sensitivity surge duration (minutes, ~360 min = 6 hours).
    pub nocturnal_surge_duration_min: f64,
}

impl Default for ExerciseParams {
    fn default() -> Self {
        Self {
            nimgu_rate_per_met: 0.22,
            glycogen_depletion_rate: 0.008,
            sensitivity_boost_max: 1.65,
            tau_acute_recovery_min: 150.0,
            nocturnal_surge_onset_min: 480.0,
            nocturnal_surge_duration_min: 360.0,
        }
    }
}

/// Dynamic physiological state of exercise exertion and post-workout adaptation.
#[derive(Debug, Clone, Default)]
pub struct ExerciseState {
    /// Current instantaneous exertion intensity (METs: 1.0 = resting).
    pub current_mets: f64,
    /// Minutes elapsed since the end of the last exercise session.
    pub time_since_exercise_min: f64,
    /// Acute post-workout insulin sensitivity multiplier (1.0 = baseline).
    pub acute_sensitivity_multiplier: f64,
    /// Delayed secondary nocturnal sensitivity multiplier.
    pub delayed_sensitivity_multiplier: f64,
    /// Total cumulative exercise expenditure (MET-minutes).
    pub cumulative_met_minutes: f64,
}

/// Exercise metabolic simulation engine.
#[derive(Debug, Clone)]
pub struct ExerciseSystem {
    pub params: ExerciseParams,
    pub state: ExerciseState,
    current_event: Option<ExerciseEvent>,
}

impl ExerciseSystem {
    /// Construct a new exercise system.
    pub fn new(params: ExerciseParams) -> Self {
        Self {
            params,
            state: ExerciseState {
                current_mets: 1.0,
                time_since_exercise_min: 9999.0,
                acute_sensitivity_multiplier: 1.0,
                delayed_sensitivity_multiplier: 1.0,
                cumulative_met_minutes: 0.0,
            },
            current_event: None,
        }
    }

    /// Schedule an exercise session.
    pub fn trigger_exercise(&mut self, event: ExerciseEvent) {
        self.current_event = Some(event);
    }

    /// Step exercise state across interval `dt` (min) at simulation time `t` (min).
    ///
    /// Returns:
    /// - Active muscle non-insulin-mediated glucose uptake (mmol/min).
    /// - Composite insulin sensitivity multiplier ($S_I$ factor $\ge 1.0$).
    pub fn step(&mut self, t: f64, dt: f64) -> (f64, f64) {
        let p = &self.params;
        let s = &mut self.state;

        let mut active_mets = 1.0;
        let mut is_exercising = false;

        if let Some(ref event) = self.current_event {
            if t >= event.start_time_min && t < event.start_time_min + event.duration_min {
                active_mets = event.intensity_mets;
                is_exercising = true;
            } else if t >= event.start_time_min + event.duration_min && is_exercising {
                s.time_since_exercise_min = 0.0;
            }
        }

        s.current_mets = active_mets;

        let nimgu_uptake = if active_mets > 1.0 {
            let delta_mets = active_mets - 1.0;
            s.cumulative_met_minutes += delta_mets * dt;
            s.time_since_exercise_min = 0.0;

            // Direct muscular uptake via mechanical contraction (GLUT4 translocation without insulin)
            delta_mets * p.nimgu_rate_per_met
        } else {
            s.time_since_exercise_min += dt;
            0.0
        };

        // Acute insulin sensitivity enhancement during and immediately post-exercise
        if is_exercising {
            let boost_target =
                1.0 + (p.sensitivity_boost_max - 1.0) * ((active_mets - 1.0) / 7.0).clamp(0.0, 1.0);
            s.acute_sensitivity_multiplier = boost_target;
        } else if s.time_since_exercise_min < p.tau_acute_recovery_min * 4.0 {
            // Exponential decay back to 1.0
            let decay = (-dt / p.tau_acute_recovery_min).exp();
            s.acute_sensitivity_multiplier = 1.0 + (s.acute_sensitivity_multiplier - 1.0) * decay;
        } else {
            s.acute_sensitivity_multiplier = 1.0;
        }

        // Delayed secondary nocturnal insulin sensitivity surge (biphasic exercise response)
        // Known clinical risk factor for severe delayed nocturnal hypoglycemia (often 8-14 hours post-workout)
        if s.cumulative_met_minutes > 60.0
            && s.time_since_exercise_min >= p.nocturnal_surge_onset_min
            && s.time_since_exercise_min
                <= p.nocturnal_surge_onset_min + p.nocturnal_surge_duration_min
        {
            let surge_progress = (s.time_since_exercise_min - p.nocturnal_surge_onset_min)
                / p.nocturnal_surge_duration_min;
            // Bell curve profile over the window
            let bell = (std::f64::consts::PI * surge_progress).sin();
            let max_delayed_boost = 0.35 * (s.cumulative_met_minutes / 200.0).clamp(0.2, 1.0);
            s.delayed_sensitivity_multiplier = 1.0 + max_delayed_boost * bell;
        } else {
            s.delayed_sensitivity_multiplier = 1.0;
        }

        let composite_sensitivity_factor =
            s.acute_sensitivity_multiplier * s.delayed_sensitivity_multiplier;

        (nimgu_uptake, composite_sensitivity_factor)
    }
}
