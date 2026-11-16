//! CLI parsing and subcommand execution unit tests.

use clap::Parser;
use glyco_soma::cli::args::{Cli, Commands};

#[test]
fn test_cli_parse_simulate_defaults() {
    let args = vec!["glyco-soma", "simulate"];
    let cli = Cli::try_parse_from(args).unwrap();

    match cli.command {
        Commands::Simulate(sim) => {
            assert_eq!(sim.duration_hours, 24.0);
            assert_eq!(sim.meals, "45,75,50");
            assert!(sim.plot);
        }
        _ => panic!("Expected Simulate subcommand"),
    }
}

#[test]
fn test_cli_parse_calibrate() {
    let args = vec![
        "glyco-soma",
        "calibrate",
        "--points",
        "40",
        "--carbs",
        "75.0",
        "--true-isf",
        "50.0",
    ];
    let cli = Cli::try_parse_from(args).unwrap();

    match cli.command {
        Commands::Calibrate(cal) => {
            assert_eq!(cal.points, 40);
            assert_eq!(cal.carbs, 75.0);
            assert_eq!(cal.true_isf, 50.0);
        }
        _ => panic!("Expected Calibrate subcommand"),
    }
}

#[test]
fn test_cli_parse_benchmark() {
    let args = vec![
        "glyco-soma",
        "benchmark",
        "--cohort-size",
        "15",
        "--days",
        "2",
        "--seed",
        "999",
    ];
    let cli = Cli::try_parse_from(args).unwrap();

    match cli.command {
        Commands::Benchmark(bm) => {
            assert_eq!(bm.cohort_size, 15);
            assert_eq!(bm.days, 2);
            assert_eq!(bm.seed, 999);
        }
        _ => panic!("Expected Benchmark subcommand"),
    }
}

#[test]
fn test_cli_parse_inspect() {
    let args = vec!["glyco-soma", "inspect", "--cohort", "pediatric"];
    let cli = Cli::try_parse_from(args).unwrap();

    match cli.command {
        Commands::Inspect(insp) => {
            assert_eq!(
                insp.cohort,
                glyco_soma::cli::args::CliPatientCohort::Pediatric
            );
        }
        _ => panic!("Expected Inspect subcommand"),
    }
}

#[test]
fn test_cli_parse_demo() {
    let args = vec!["glyco-soma", "demo", "--cohort", "adolescent"];
    let cli = Cli::try_parse_from(args).unwrap();

    match cli.command {
        Commands::Demo(demo) => {
            assert_eq!(
                demo.cohort,
                glyco_soma::cli::args::CliPatientCohort::Adolescent
            );
        }
        _ => panic!("Expected Demo subcommand"),
    }
}
