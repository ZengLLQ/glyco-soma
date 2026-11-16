//! Continuous Glucose Monitor (CGM) sensors, noise artifacts, and Clarke Error Grid clinical accuracy validation.

pub mod cgm;
pub mod clark_grid;

pub use cgm::{CgmParams, CgmSensor, CgmState};
pub use clark_grid::{classify_clark_point, evaluate_clark_grid, ClarkGridReport, ClarkZone};
