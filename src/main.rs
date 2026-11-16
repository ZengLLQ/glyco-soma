//! CLI binary entrypoint for GlycoSoma.

#![allow(clippy::approx_constant, clippy::manual_str_repeat)]

use clap::Parser;
use glyco_soma::cli::args::{Cli, Commands};
use glyco_soma::cli::runner::{run_benchmark, run_calibrate, run_demo, run_inspect, run_simulate};
use glyco_soma::error::Result;

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Simulate(args) => run_simulate(args)?,
        Commands::Calibrate(args) => run_calibrate(args)?,
        Commands::Benchmark(args) => run_benchmark(args)?,
        Commands::Inspect(args) => run_inspect(args)?,
        Commands::Demo(args) => run_demo(args)?,
    }

    Ok(())
}
