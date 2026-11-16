//! # GlycoSoma
//!
//! **Biophysical simulation engine for multi-hormonal glucose-insulin-glucagon homeostasis, delayed subcutaneous absorption kinetics, and adaptive closed-loop artificial pancreas control.**
//!
//! ## Subsystems
//! - **Mathematical Engines (`math`)**: Adaptive Runge-Kutta-Fehlberg 4(5) (RKF45), Crank-Nicolson stiff implicit solver, and delay differential equation (DDE) Hermite cubic ring buffer.
//! - **Endocrine & Digestive Models (`models`)**: Hovorka 2004 multi-compartment non-linear model, Bergman Minimal Model, 2-compartment gastric emptying with intestinal absorption feedback, counter-regulatory glucagon kinetics, and exercise-induced metabolic exertion (NIMGU).
//! - **Sensors & Error Grids (`sensors`)**: Continuous Glucose Monitor (CGM) interstitial lag, colored AR(1) noise, compression artifacts, and Clarke Error Grid Analysis (CEGA) with exact polygon point-in-polygon ray casting.
//! - **Control Systems (`controllers`)**: Constrained Model Predictive Control (MPC) with asymmetric hypoglycemia barrier penalties, dual-hormone PID with anti-windup clamping, Predictive Low Glucose Suspend (PLGS) state machines, and traditional open-loop basal-bolus therapy.
//! - **Clinical Metrics & Virtual Cohorts (`clinical`)**: International ATTD 2024 consensus reports (TIR, TBR, TAR, GMI, GRI), glycemic variability (MAGE, CV%, LBGI, HBGI), and 100-patient in-silico cohort generation.
//! - **7-Bit ASCII Visualizers (`visualizer`)**: Strip chart plotting, 24-hour Ambulatory Glucose Profile (AGP) modal reports, and Clarke Error Grid scatter diagrams using strictly ASCII codes 0–127.

#![allow(
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::approx_constant,
    clippy::useless_vec,
    clippy::manual_str_repeat
)]

pub mod cli;
pub mod clinical;
pub mod controllers;
pub mod error;
pub mod math;
pub mod models;
pub mod sensors;
pub mod types;
pub mod visualizer;

pub mod prelude {
    pub use crate::clinical::{
        compute_attd_metrics, compute_variability_metrics, AttdMetricsReport,
        GlycemicVariabilityReport, VirtualCohortGenerator,
    };
    pub use crate::controllers::{
        BasalBolusController, ModelPredictiveController, MpcConfig, PidConfig, PidController,
        PlgsConfig, PlgsController,
    };
    pub use crate::error::{GlycoError, Result};
    pub use crate::math::{Rkf45Config, Rkf45Solver, StiffSolver, StiffSolverConfig};
    pub use crate::models::{
        CounterRegulationSystem, DigestionParams, DigestionSystem, ExerciseParams, ExerciseSystem,
        HovorkaParams, HovorkaSystem, IntegratedMetabolicEngine,
    };
    pub use crate::sensors::{
        classify_clark_point, evaluate_clark_grid, CgmParams, CgmSensor, ClarkGridReport, ClarkZone,
    };
    pub use crate::types::{
        ExerciseEvent, HormoneInfusion, MealEvent, MetabolicState, PatientCohort, PatientProfile,
    };
    pub use crate::visualizer::{
        render_ascii_agp_report, render_ascii_clark_plot, render_ascii_glucose_chart,
    };
}
