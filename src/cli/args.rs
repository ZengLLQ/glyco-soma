//! Command-line argument definitions and subcommands for GlycoSoma.

use clap::{Args, Parser, Subcommand, ValueEnum};

/// Top-level CLI argument structure.
#[derive(Parser, Debug)]
#[command(
    name = "glyco-soma",
    about = "Biophysical simulation engine for multi-hormonal glucose-insulin-glucagon homeostasis and closed-loop artificial pancreas control.",
    version = "0.1.0"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

/// Available operational subcommands.
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Run an in-silico virtual patient simulation across meals, exercise, and closed-loop control.
    Simulate(SimulateArgs),
    /// Calibrate patient-specific insulin sensitivity and carb ratio from synthetic CGM traces.
    Calibrate(CalibrateArgs),
    /// Run a comparative head-to-head clinical benchmark comparing control modalities across a cohort.
    Benchmark(BenchmarkArgs),
    /// Inspect virtual patient metabolic characteristics, default profiles, and model equations.
    Inspect(InspectArgs),
    /// Run a rich live interactive demonstration showing ASCII glucose trajectories and ATTD reports.
    Demo(DemoArgs),
}

/// Controller type selection for CLI.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum CliControllerMode {
    /// Traditional open-loop basal-bolus pump therapy.
    OpenLoop,
    /// Automated Predictive Low Glucose Suspend (PLGS).
    Plgs,
    /// Dual-hormone feedback PID controller.
    Pid,
    /// Constrained Model Predictive Control (MPC) with asymmetric cost functions.
    Mpc,
}

/// Virtual patient cohort selection for CLI.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum CliPatientCohort {
    /// Adult patient profile (70 kg).
    Adult,
    /// Adolescent profile with pubertal resistance (55 kg).
    Adolescent,
    /// Pediatric profile with high sensitivity (28 kg).
    Pediatric,
}

/// Arguments for the `simulate` subcommand.
#[derive(Args, Debug)]
pub struct SimulateArgs {
    /// Patient cohort selection.
    #[arg(short = 'p', long, value_enum, default_value = "adult")]
    pub cohort: CliPatientCohort,

    /// Controller modality selection.
    #[arg(short = 'c', long, value_enum, default_value = "mpc")]
    pub controller: CliControllerMode,

    /// Total simulation duration in hours (e.g. 24.0 or 48.0).
    #[arg(short = 'd', long, default_value_t = 24.0)]
    pub duration_hours: f64,

    /// Meal carb amounts in grams (e.g. 45,70,50).
    #[arg(long, default_value = "45,75,50")]
    pub meals: String,

    /// Meal times in hours (e.g. 7.5,12.5,18.5).
    #[arg(long, default_value = "7.5,12.5,18.5")]
    pub meal_times: String,

    /// Add an exercise session (onset hour and duration in minutes, e.g. 15.0:45).
    #[arg(long)]
    pub exercise: Option<String>,

    /// Save full trajectory snapshot to JSON file.
    #[arg(short, long)]
    pub output: Option<String>,

    /// Render 7-bit ASCII glucose strip chart in terminal.
    #[arg(long, default_value_t = true)]
    pub plot: bool,
}

/// Arguments for the `calibrate` subcommand.
#[derive(Args, Debug)]
pub struct CalibrateArgs {
    /// Number of synthetic post-prandial CGM points to generate and fit.
    #[arg(short, long, default_value_t = 36)]
    pub points: usize,

    /// Meal size in grams of carbohydrates used for calibration trial.
    #[arg(long, default_value_t = 60.0)]
    pub carbs: f64,

    /// True patient ISF (mg/dL/U) for synthetic generation.
    #[arg(long, default_value_t = 45.0)]
    pub true_isf: f64,

    /// True patient CIR (g/U) for synthetic generation.
    #[arg(long, default_value_t = 10.0)]
    pub true_cir: f64,
}

/// Arguments for the `benchmark` subcommand.
#[derive(Args, Debug)]
pub struct BenchmarkArgs {
    /// Number of virtual patients in trial cohort (e.g. 20 to 100).
    #[arg(short, long, default_value_t = 25)]
    pub cohort_size: usize,

    /// Trial duration in days per patient (e.g. 1 to 7 days).
    #[arg(short, long, default_value_t = 1)]
    pub days: usize,

    /// Random seed for trial reproducibility.
    #[arg(long, default_value_t = 1337)]
    pub seed: u64,
}

/// Arguments for the `inspect` subcommand.
#[derive(Args, Debug)]
pub struct InspectArgs {
    /// Patient cohort to inspect.
    #[arg(short, long, value_enum, default_value = "adult")]
    pub cohort: CliPatientCohort,
}

/// Arguments for the `demo` subcommand.
#[derive(Args, Debug)]
pub struct DemoArgs {
    /// Patient cohort for demo.
    #[arg(short, long, value_enum, default_value = "adult")]
    pub cohort: CliPatientCohort,
}
