//! Mathematical foundations, numerical ODE solvers, delay differential equations, and statistics.

pub mod delay_buffer;
pub mod ode_rkf45;
pub mod ode_stiff;
pub mod statistics;

pub use delay_buffer::DelayBuffer;
pub use ode_rkf45::{Rkf45Config, Rkf45Solver, SolverStats};
pub use ode_stiff::{StiffSolver, StiffSolverConfig};
pub use statistics::{
    compute_mard, compute_pearson_r, compute_rmse, compute_series_stats, percentile_sorted,
    GaussianKde, SeriesStats,
};
