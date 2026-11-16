//! Integrated endocrine, digestive, counter-regulatory, and exercise metabolic simulation models.

pub mod bergman;
pub mod counter_regulation;
pub mod digestion;
pub mod exercise;
pub mod hovorka;

pub use bergman::{BergmanParams, BergmanSystem};
pub use counter_regulation::{CounterRegulationSystem, GlucagonParams, GlucagonState};
pub use digestion::{DigestionParams, DigestionState, DigestionSystem};
pub use exercise::{ExerciseParams, ExerciseState, ExerciseSystem};
pub use hovorka::{state_idx, HovorkaParams, HovorkaSystem, HOVORKA_STATE_DIM};

use crate::error::Result;
use crate::math::ode_rkf45::Rkf45Solver;
use crate::types::{HormoneInfusion, MetabolicState, PatientProfile};

/// Integrated whole-body physiological simulation model coupling Hovorka endocrine kinetics,
/// digestion/gastric emptying, counter-regulatory glucagon dynamics, and physical exercise.
#[derive(Debug, Clone)]
pub struct IntegratedMetabolicEngine {
    pub patient: PatientProfile,
    pub hovorka: HovorkaSystem,
    pub digestion: DigestionSystem,
    pub counter_reg: CounterRegulationSystem,
    pub exercise: ExerciseSystem,
    pub state_vector: Vec<f64>,
}

impl IntegratedMetabolicEngine {
    /// Construct a new integrated whole-body metabolic engine from a patient profile.
    pub fn new(patient: PatientProfile) -> Self {
        let hovorka_params = HovorkaParams::from_patient(&patient);
        let hovorka = HovorkaSystem::new(hovorka_params);
        let digestion = DigestionSystem::new(DigestionParams::default());
        let counter_reg = CounterRegulationSystem::new(GlucagonParams::default());
        let exercise = ExerciseSystem::new(ExerciseParams::default());

        let initial_state = MetabolicState::initial_basal(&patient);

        // State vector matching HOVORKA_STATE_DIM
        let mut state_vector = vec![0.0; HOVORKA_STATE_DIM];
        state_vector[state_idx::S1] = initial_state.s1_subcut_insulin;
        state_vector[state_idx::S2] = initial_state.s2_subcut_insulin;
        state_vector[state_idx::I] = initial_state.insulin_plasma_mu_l * 6.0; // pmol/L
        state_vector[state_idx::X1] = initial_state.x1_remote_transport;
        state_vector[state_idx::X2] = initial_state.x2_remote_disposal;
        state_vector[state_idx::X3] = initial_state.x3_remote_egp;
        state_vector[state_idx::Q1] = initial_state.q1_accessible_glucose_mmol;
        state_vector[state_idx::Q2] = initial_state.q2_nonaccessible_glucose_mmol;

        Self {
            patient,
            hovorka,
            digestion,
            counter_reg,
            exercise,
            state_vector,
        }
    }

    /// Advance the complete coupled whole-body system across time step `dt` (min) at simulation time `t` (min).
    pub fn step(
        &mut self,
        t: f64,
        dt: f64,
        infusion: &HormoneInfusion,
        solver: &mut Rkf45Solver,
    ) -> Result<MetabolicState> {
        // 1. Advance digestive absorption
        let gut_absorption_rate_mmol_min = self.digestion.step(dt);

        // 2. Advance exercise exertion
        let (exercise_uptake_mmol_min, exercise_sensitivity_multiplier) = self.exercise.step(t, dt);

        // Modulate Hovorka insulin sensitivity dynamically based on exercise
        let current_plasma_glucose = self.hovorka.plasma_glucose_mg_dl(&self.state_vector);
        let current_plasma_insulin = self.hovorka.plasma_insulin_mu_l(&self.state_vector);

        // 3. Advance glucagon counter-regulation
        let glucagon_egp_boost_mmol_min = self.counter_reg.step(
            dt,
            infusion,
            self.patient.body_weight_kg,
            current_plasma_glucose,
            current_plasma_insulin,
        );

        // Apply exercise sensitivity multiplier to effective sensitivity parameters
        let base_params = HovorkaParams::from_patient(&self.patient);
        self.hovorka.params.k_b1 = base_params.k_b1 * exercise_sensitivity_multiplier;
        self.hovorka.params.k_b2 = base_params.k_b2 * exercise_sensitivity_multiplier;

        // Dawn phenomenon circadian modulation on endogenous glucose production (04:00-08:00 = 240-480 min)
        let minute_of_day = (t % 1440.0).max(0.0);
        if (240.0..=480.0).contains(&minute_of_day) {
            let dawn_phase = ((minute_of_day - 240.0) / 240.0) * std::f64::consts::PI;
            let dawn_surge = dawn_phase.sin() * self.patient.dawn_phenomenon_amp;
            self.hovorka.params.egp_0 = base_params.egp_0 * (1.0 + dawn_surge);
        } else {
            self.hovorka.params.egp_0 = base_params.egp_0;
        }

        // 4. Deliver discrete insulin bolus directly into subcutaneous depot S1
        self.state_vector[state_idx::S1] +=
            infusion.insulin_bolus_u * crate::types::INSULIN_U_TO_PMOL;

        // 5. Integrate Hovorka differential equations using adaptive RKF45
        let hovorka_ref = &self.hovorka;
        let infusion_copy = *infusion;

        solver.integrate_step(t, t + dt, &mut self.state_vector, |_time, y, dydt| {
            let _ = hovorka_ref.evaluate_derivatives(
                _time,
                y,
                &infusion_copy,
                gut_absorption_rate_mmol_min,
                exercise_uptake_mmol_min,
                glucagon_egp_boost_mmol_min,
                dydt,
            );
        })?;

        // 5. Package full metabolic state snapshot
        let plasma_g = self.hovorka.plasma_glucose_mg_dl(&self.state_vector);
        let plasma_i = self.hovorka.plasma_insulin_mu_l(&self.state_vector);

        Ok(MetabolicState {
            time_min: t + dt,
            glucose_plasma_mg_dl: plasma_g,
            glucose_interstitial_mg_dl: plasma_g, // updated by sensor module
            insulin_plasma_mu_l: plasma_i,
            s1_subcut_insulin: self.state_vector[state_idx::S1],
            s2_subcut_insulin: self.state_vector[state_idx::S2],
            x1_remote_transport: self.state_vector[state_idx::X1],
            x2_remote_disposal: self.state_vector[state_idx::X2],
            x3_remote_egp: self.state_vector[state_idx::X3],
            q1_accessible_glucose_mmol: self.state_vector[state_idx::Q1],
            q2_nonaccessible_glucose_mmol: self.state_vector[state_idx::Q2],
            stomach_solid_g: self.digestion.state.stomach_solid_g,
            stomach_liquid_g: self.digestion.state.stomach_liquid_g,
            intestine_g: self.digestion.state.gut_proximal_g + self.digestion.state.gut_distal_g,
            glucagon_plasma_pg_ml: self.counter_reg.state.plasma_glucagon_pg_ml,
            liver_glycogen_fraction: self.counter_reg.state.glycogen_fraction,
            active_exercise_uptake_mmol_min: exercise_uptake_mmol_min,
        })
    }
}
