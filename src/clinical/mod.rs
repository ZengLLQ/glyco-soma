//! Clinical standards, ATTD consensus metrics, glycemic variability, and in-silico population trials.

pub mod attd_metrics;
pub mod calibration;
pub mod variability;
pub mod virtual_cohort;

pub use attd_metrics::{compute_attd_metrics, AttdMetricsReport};
pub use calibration::{calibrate_patient_parameters, CalibrationResult};
pub use variability::{compute_variability_metrics, GlycemicVariabilityReport};
pub use virtual_cohort::VirtualCohortGenerator;
