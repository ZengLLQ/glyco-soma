//! Error types and result definitions for GlycoSoma.

use std::fmt;

/// Primary error type for GlycoSoma operations.
#[derive(Debug)]
pub enum GlycoError {
    /// Provided numerical parameter is out of physiological bounds.
    InvalidParameter {
        parameter: String,
        value: f64,
        reason: String,
    },
    /// Numerical ODE solver divergence or failure to converge.
    SolverDivergence {
        step: usize,
        time: f64,
        error_norm: f64,
    },
    /// Negative biological mass or concentration detected where strictly positive is required.
    NegativeMassDetected {
        compartment: String,
        value: f64,
        time: f64,
    },
    /// Time step reduced below numerical stability threshold.
    TimeStepTooSmall { dt: f64, min_dt: f64, time: f64 },
    /// Matrix inversion or linear algebra failure during controller optimization.
    MatrixSingular { dimension: usize, detail: String },
    /// Optimization or parameter fitting failed to find optimal solution.
    OptimizationFailed {
        iterations: usize,
        residual: f64,
        reason: String,
    },
    /// Empty data series provided for statistical or clinical analysis.
    EmptyDataSeries { series_name: String },
    /// Standard I/O error.
    Io(std::io::Error),
    /// JSON serialization/deserialization error.
    Json(serde_json::Error),
}

impl fmt::Display for GlycoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GlycoError::InvalidParameter {
                parameter,
                value,
                reason,
            } => {
                write!(
                    f,
                    "Invalid parameter '{}' = {}: {}",
                    parameter, value, reason
                )
            }
            GlycoError::SolverDivergence {
                step,
                time,
                error_norm,
            } => {
                write!(
                    f,
                    "Solver diverged at step {} (t = {:.2} min) with error norm {:.2e}",
                    step, time, error_norm
                )
            }
            GlycoError::NegativeMassDetected {
                compartment,
                value,
                time,
            } => {
                write!(
                    f,
                    "Negative mass in compartment '{}' = {:.4} at t = {:.2} min",
                    compartment, value, time
                )
            }
            GlycoError::TimeStepTooSmall { dt, min_dt, time } => {
                write!(
                    f,
                    "Adaptive time step dt = {:.2e} fell below minimum {:.2e} at t = {:.2} min",
                    dt, min_dt, time
                )
            }
            GlycoError::MatrixSingular { dimension, detail } => {
                write!(
                    f,
                    "Matrix of dimension {}x{} is singular or non-invertible: {}",
                    dimension, dimension, detail
                )
            }
            GlycoError::OptimizationFailed {
                iterations,
                residual,
                reason,
            } => {
                write!(
                    f,
                    "Optimization failed after {} iterations (residual = {:.4e}): {}",
                    iterations, residual, reason
                )
            }
            GlycoError::EmptyDataSeries { series_name } => {
                write!(
                    f,
                    "Data series '{}' is empty, cannot compute statistics",
                    series_name
                )
            }
            GlycoError::Io(err) => write!(f, "I/O error: {}", err),
            GlycoError::Json(err) => write!(f, "JSON error: {}", err),
        }
    }
}

impl std::error::Error for GlycoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            GlycoError::Io(err) => Some(err),
            GlycoError::Json(err) => Some(err),
            _ => None,
        }
    }
}

impl From<std::io::Error> for GlycoError {
    fn from(err: std::io::Error) -> Self {
        GlycoError::Io(err)
    }
}

impl From<serde_json::Error> for GlycoError {
    fn from(err: serde_json::Error) -> Self {
        GlycoError::Json(err)
    }
}

/// Convenience result alias for GlycoSoma operations.
pub type Result<T> = std::result::Result<T, GlycoError>;
