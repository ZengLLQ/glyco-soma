//! Command-line interface definitions and command runners.

pub mod args;
pub mod runner;

pub use args::{Cli, Commands};
pub use runner::{run_benchmark, run_calibrate, run_demo, run_inspect, run_simulate};
