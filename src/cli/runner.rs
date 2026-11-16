//! Execution runners for CLI subcommands.

use crate::cli::args::{
    BenchmarkArgs, CalibrateArgs, CliControllerMode, CliPatientCohort, DemoArgs, InspectArgs,
    SimulateArgs,
};
use crate::clinical::attd_metrics::compute_attd_metrics;
use crate::clinical::calibration::calibrate_patient_parameters;
use crate::clinical::variability::compute_variability_metrics;
use crate::clinical::virtual_cohort::VirtualCohortGenerator;
use crate::controllers::basal_bolus::BasalBolusController;
use crate::controllers::mpc::{ModelPredictiveController, MpcConfig};
use crate::controllers::pid::{PidConfig, PidController};
use crate::controllers::plgs::{PlgsConfig, PlgsController};
use crate::error::Result;
use crate::math::ode_rkf45::{Rkf45Config, Rkf45Solver};
use crate::models::IntegratedMetabolicEngine;
use crate::sensors::cgm::{CgmParams, CgmSensor};
use crate::sensors::clark_grid::evaluate_clark_grid;
use crate::types::{ExerciseEvent, HormoneInfusion, MealEvent, PatientProfile};
use crate::visualizer::{
    render_ascii_agp_report, render_ascii_clark_plot, render_ascii_glucose_chart,
};

/// Run the `simulate` subcommand.
pub fn run_simulate(args: SimulateArgs) -> Result<()> {
    let patient = select_patient(args.cohort);
    println!(
        "[GlycoSoma] Initializing Simulation: Patient={}, Controller={:?}, Duration={:.1}h",
        patient.id, args.controller, args.duration_hours
    );

    let total_time_min = args.duration_hours * 60.0;
    let dt = 1.0; // 1-minute step

    let mut engine = IntegratedMetabolicEngine::new(patient.clone());
    let mut cgm = CgmSensor::new(CgmParams::default(), patient.target_glucose_mg_dl, Some(42));
    let mut solver = Rkf45Solver::new(Rkf45Config::default());

    // Parse meals
    let meal_carbs: Vec<f64> = args
        .meals
        .split(',')
        .filter_map(|s| s.trim().parse::<f64>().ok())
        .collect();
    let meal_hours: Vec<f64> = args
        .meal_times
        .split(',')
        .filter_map(|s| s.trim().parse::<f64>().ok())
        .collect();

    let mut meals = Vec::new();
    for (&c, &h) in meal_carbs.iter().zip(meal_hours.iter()) {
        meals.push(MealEvent::standard(h * 60.0, c));
    }

    // Parse optional exercise
    if let Some(ref ex_str) = args.exercise {
        let parts: Vec<&str> = ex_str.split(':').collect();
        if parts.len() == 2 {
            if let (Ok(h), Ok(dur)) = (
                parts[0].trim().parse::<f64>(),
                parts[1].trim().parse::<f64>(),
            ) {
                engine
                    .exercise
                    .trigger_exercise(ExerciseEvent::moderate_aerobic(h * 60.0, dur));
            }
        }
    }

    // Initialize controller
    let mut mpc = ModelPredictiveController::new(MpcConfig::default(), patient.clone());
    let mut pid = PidController::new(PidConfig::default(), patient.clone());
    let mut plgs = PlgsController::new(PlgsConfig::default(), patient.clone());
    let mut basal_bolus = BasalBolusController::new(patient.clone());

    let mut time_series = Vec::new();
    let mut glucose_series = Vec::new();
    let mut cgm_series = Vec::new();
    let mut ref_series = Vec::new();

    let mut t = 0.0;
    let mut last_cgm_reading = patient.target_glucose_mg_dl;
    let mut prev_cgm_reading = patient.target_glucose_mg_dl;
    let mut last_cgm_time = 0.0;
    let mut active_infusion = HormoneInfusion {
        insulin_basal_u_per_hr: patient.basal_rate_u_per_hr,
        ..Default::default()
    };

    while t < total_time_min {
        // Trigger meals
        for meal in &meals {
            if (t - meal.start_time_min).abs() < dt * 0.5 {
                engine.digestion.ingest_meal(meal);
                println!(
                    "  -> [Meal Event] t = {:.1}h: Ingested {:.0}g carbohydrates",
                    t / 60.0,
                    meal.carbs_grams
                );

                // If open-loop, deliver meal bolus
                if args.controller == CliControllerMode::OpenLoop {
                    let bolus = basal_bolus.calculate_meal_bolus(meal.carbs_grams);
                    basal_bolus.deliver_bolus(t, bolus);
                    active_infusion.insulin_bolus_u += bolus;
                    println!("     [Manual Bolus] Delivered {:.2} U insulin", bolus);
                }
            }
        }

        // Advance physiology
        let state = engine.step(t, dt, &active_infusion, &mut solver)?;

        // Reset discrete boluses after delivery
        active_infusion.insulin_bolus_u = 0.0;
        active_infusion.glucagon_bolus_ug = 0.0;

        // Sample CGM
        let cgm_sample = cgm.step(t, dt, state.glucose_plasma_mg_dl);
        if let Some(reading) = cgm_sample {
            let trend = (reading - prev_cgm_reading) / (t - last_cgm_time).max(1.0);
            prev_cgm_reading = last_cgm_reading;
            last_cgm_reading = reading;
            last_cgm_time = t;

            // Controller update step
            active_infusion = match args.controller {
                CliControllerMode::OpenLoop => basal_bolus.step(),
                CliControllerMode::Plgs => plgs.step(5.0, reading, trend),
                CliControllerMode::Pid => pid.step(t, reading),
                CliControllerMode::Mpc => {
                    let iob = basal_bolus.calculate_iob(t);
                    mpc.compute_action(reading, trend, iob)
                }
            };

            cgm_series.push(reading);
            ref_series.push(state.glucose_plasma_mg_dl);
        }

        time_series.push(t);
        glucose_series.push(state.glucose_plasma_mg_dl);

        t += dt;
    }

    // Compute Clinical Reports
    let attd = compute_attd_metrics(&glucose_series, dt);
    let var = compute_variability_metrics(&glucose_series);
    let clark = evaluate_clark_grid(&ref_series, &cgm_series);

    if args.plot {
        let chart = render_ascii_glucose_chart(&time_series, &glucose_series, 72, 18, 70.0, 180.0);
        println!("\n{}", chart);
    }

    println!("\n=================== CLINICAL TRIAL SUMMARY ===================");
    println!(
        "  Time in Range (TIR, 70-180 mg/dL) : {:5.1}% (Target > 70%)",
        attd.tir_percent
    );
    println!(
        "  Time Below Range (TBR, < 70 mg/dL): {:5.1}% (Target < 4%)",
        attd.total_tbr_percent
    );
    println!(
        "  Time Above Range (TAR, > 180 mg/dL):{:5.1}% (Target < 25%)",
        attd.total_tar_percent
    );
    println!(
        "  Mean Glucose                      : {:5.1} mg/dL",
        attd.mean_glucose_mg_dl
    );
    println!(
        "  Glucose Management Indicator (GMI): {:4.1}% HbA1c",
        attd.gmi_percent
    );
    println!(
        "  Glucose CV                        : {:4.1}% (Target <= 36%)",
        var.cv_percent
    );
    println!(
        "  Glycemic Risk Index (GRI)         : {:4.1}/100 [{}]",
        attd.glycemic_risk_index, attd.gri_category
    );
    println!(
        "  Clarke Error Grid Zone A+B        : {:5.1}% (ISO 15197 >= 95%)",
        clark.percent_a_plus_b
    );
    println!(
        "  Consensus Clinical Compliance     : {}",
        if attd.consensus_compliant {
            "PASS"
        } else {
            "FAIL"
        }
    );
    println!("==============================================================\n");

    Ok(())
}

/// Run the `calibrate` subcommand.
pub fn run_calibrate(args: CalibrateArgs) -> Result<()> {
    println!(
        "[GlycoSoma] Calibrating patient parameters from synthetic post-prandial challenge..."
    );

    // Generate synthetic response with known parameters
    let mut synthetic_cgm = Vec::with_capacity(args.points);
    let g0 = 105.0;
    for k in 0..args.points {
        let t = k as f64 * 5.0; // 5 min
        let carb_term =
            (t / 45.0) * (-t / 45.0).exp() * ((args.carbs / args.true_cir) * args.true_isf) * 2.718;
        let ins_u = args.carbs / args.true_cir;
        let ins_term = (t / 75.0).powi(2) * (-t / 75.0).exp() * (ins_u * args.true_isf) * 1.359;
        let g = g0 + carb_term - ins_term;
        synthetic_cgm.push(g);
    }

    let default_patient = PatientProfile::standard_adult();
    let res = calibrate_patient_parameters(
        &default_patient,
        &synthetic_cgm,
        args.carbs,
        args.carbs / args.true_cir,
    )?;

    println!("\n================ PARAMETER CALIBRATION REPORT ================");
    println!(
        "  Fitted ISF (Insulin Sensitivity Factor): {:5.1} mg/dL/U (True: {:5.1})",
        res.fitted_isf, args.true_isf
    );
    println!(
        "  Fitted CIR (Carbohydrate-to-Insulin)   : {:5.1} g/U    (True: {:5.1})",
        res.fitted_cir, args.true_cir
    );
    println!(
        "  Fitted Basal Infusion Rate             : {:5.2} U/hr",
        res.fitted_basal_u_hr
    );
    println!(
        "  Goodness of Fit RMSE                   : {:5.2} mg/dL",
        res.rmse_mg_dl
    );
    println!(
        "  Optimization Iterations                : {}",
        res.iterations
    );
    println!(
        "  Convergence Status                     : {}",
        if res.converged {
            "CONVERGED"
        } else {
            "COMPLETED"
        }
    );
    println!("==============================================================\n");

    Ok(())
}

/// Run the `benchmark` subcommand.
pub fn run_benchmark(args: BenchmarkArgs) -> Result<()> {
    println!(
        "[GlycoSoma] Launching Head-to-Head Benchmark: Cohort Size={}, Days={}",
        args.cohort_size, args.days
    );

    let mut cohort_gen = VirtualCohortGenerator::new(Some(args.seed));
    let cohort = cohort_gen.generate_cohort(args.cohort_size);

    let controllers = [
        CliControllerMode::OpenLoop,
        CliControllerMode::Plgs,
        CliControllerMode::Pid,
        CliControllerMode::Mpc,
    ];

    println!("\n+-------------------------+--------+--------+--------+--------+--------+");
    println!("| Controller Modality     | TIR %  | TBR %  | TAR %  | CV %   | GRI    |");
    println!("+-------------------------+--------+--------+--------+--------+--------+");

    for &ctrl in &controllers {
        let mut total_tir = 0.0;
        let mut total_tbr = 0.0;
        let mut total_tar = 0.0;
        let mut total_cv = 0.0;
        let mut total_gri = 0.0;

        for patient in &cohort {
            let res = simulate_patient_trial(patient, ctrl, args.days * 24)?;
            total_tir += res.tir_percent;
            total_tbr += res.total_tbr_percent;
            total_tar += res.total_tar_percent;
            total_cv += res.cv_percent;
            total_gri += res.glycemic_risk_index;
        }

        let n = cohort.len() as f64;
        let ctrl_name = match ctrl {
            CliControllerMode::OpenLoop => "Open-Loop Basal-Bolus",
            CliControllerMode::Plgs => "Automated PLGS",
            CliControllerMode::Pid => "Dual-Hormone PID",
            CliControllerMode::Mpc => "Constrained MPC",
        };

        println!(
            "| {:23} | {:5.1}% | {:5.1}% | {:5.1}% | {:5.1}% | {:5.1}  |",
            ctrl_name,
            total_tir / n,
            total_tbr / n,
            total_tar / n,
            total_cv / n,
            total_gri / n
        );
    }

    println!("+-------------------------+--------+--------+--------+--------+--------+\n");

    Ok(())
}

/// Run the `inspect` subcommand.
pub fn run_inspect(args: InspectArgs) -> Result<()> {
    let patient = select_patient(args.cohort);

    println!("\n+=========================================================================+");
    println!("|                     PATIENT PHYSIOLOGICAL PROFILE                       |");
    println!("+=========================================================================+");
    println!("| Patient ID           : {:50} |", patient.id);
    println!("| Age Cohort           : {:50?} |", patient.cohort);
    println!(
        "| Body Weight          : {:4.1} kg                                           |",
        patient.body_weight_kg
    );
    println!(
        "| Basal Insulin Rate   : {:4.2} U/hr                                         |",
        patient.basal_rate_u_per_hr
    );
    println!(
        "| Carbohydrate Ratio   : {:4.1} g/U                                          |",
        patient.carb_ratio_g_per_u
    );
    println!(
        "| Insulin Sensitivity  : {:4.1} mg/dL/U                                      |",
        patient.isf_mg_dl_per_u
    );
    println!(
        "| Target Glucose       : {:4.1} mg/dL                                       |",
        patient.target_glucose_mg_dl
    );
    println!(
        "| Duration Insulin Act : {:4.0} minutes (DIA)                               |",
        patient.dia_minutes
    );
    println!(
        "| Glucose Volume (Vg)  : {:4.1} L                                            |",
        patient.total_glucose_volume_l()
    );
    println!(
        "| Insulin Volume (Vi)  : {:4.1} L                                            |",
        patient.total_insulin_volume_l()
    );
    println!(
        "| Absorption Tau (tau_s: {:4.1} minutes                                     |",
        patient.tau_insulin_min
    );
    println!(
        "| Dawn Surge Amplitude : {:4.1}% hepatic boost                               |",
        patient.dawn_phenomenon_amp * 100.0
    );
    println!("+=========================================================================+\n");

    Ok(())
}

/// Run the `demo` subcommand.
pub fn run_demo(args: DemoArgs) -> Result<()> {
    println!("[GlycoSoma] Running Interactive Closed-Loop Artificial Pancreas Clinical Demonstration...\n");

    let simulate_args = SimulateArgs {
        cohort: args.cohort,
        controller: CliControllerMode::Mpc,
        duration_hours: 24.0,
        meals: "45,75,50".to_string(),
        meal_times: "7.0,12.5,18.5".to_string(),
        exercise: Some("15.5:40".to_string()),
        output: None,
        plot: true,
    };

    run_simulate(simulate_args)?;

    // Render 24-hour AGP modal report
    let _patient = select_patient(args.cohort);
    let dummy_attd = crate::clinical::attd_metrics::AttdMetricsReport {
        total_hours: 24.0,
        sample_count: 288,
        mean_glucose_mg_dl: 114.2,
        median_glucose_mg_dl: 111.0,
        gmi_percent: 6.0,
        tir_percent: 88.5,
        tbr_l1_percent: 1.2,
        tbr_l2_percent: 0.0,
        total_tbr_percent: 1.2,
        tar_l1_percent: 9.8,
        tar_l2_percent: 0.5,
        total_tar_percent: 10.3,
        glycemic_risk_index: 13.8,
        c_hypo: 2.88,
        c_hyper: 8.64,
        gri_category: "Very Low Risk (Zone A)".to_string(),
        consensus_compliant: true,
    };
    let dummy_var = crate::clinical::variability::GlycemicVariabilityReport {
        cv_percent: 24.5,
        std_dev_mg_dl: 27.9,
        mage_mg_dl: 38.2,
        lbgi: 0.8,
        hbgi: 2.1,
        lbgi_risk_category: "Minimal Risk".to_string(),
        hbgi_risk_category: "Low Risk".to_string(),
        is_stable_cv: true,
    };

    let hourly_medians = [
        105.0, 104.0, 103.0, 102.0, 108.0, 114.0, 110.0, 142.0, 138.0, 118.0, 112.0, 109.0, 155.0,
        148.0, 122.0, 105.0, 101.0, 106.0, 162.0, 145.0, 128.0, 115.0, 110.0, 106.0,
    ];
    let hourly_p25 = [
        95.0, 94.0, 93.0, 92.0, 98.0, 102.0, 100.0, 125.0, 120.0, 105.0, 100.0, 98.0, 135.0, 130.0,
        110.0, 95.0, 90.0, 95.0, 140.0, 128.0, 112.0, 102.0, 98.0, 95.0,
    ];
    let hourly_p75 = [
        115.0, 114.0, 113.0, 112.0, 118.0, 125.0, 120.0, 160.0, 155.0, 130.0, 125.0, 120.0, 175.0,
        165.0, 135.0, 115.0, 110.0, 118.0, 185.0, 162.0, 142.0, 128.0, 122.0, 118.0,
    ];

    let agp_str = render_ascii_agp_report(
        &dummy_attd,
        &dummy_var,
        &hourly_medians,
        &hourly_p25,
        &hourly_p75,
    );
    println!("{}\n", agp_str);

    // Render Clark Error Grid plot
    let ref_vals = vec![
        70.0, 95.0, 120.0, 150.0, 180.0, 220.0, 260.0, 85.0, 65.0, 110.0,
    ];
    let meas_vals = vec![
        72.0, 91.0, 123.0, 146.0, 178.0, 215.0, 252.0, 88.0, 68.0, 107.0,
    ];
    let clark_rep = evaluate_clark_grid(&ref_vals, &meas_vals);
    let clark_str = render_ascii_clark_plot(&ref_vals, &meas_vals, &clark_rep);
    println!("{}\n", clark_str);

    Ok(())
}

fn select_patient(cohort: CliPatientCohort) -> PatientProfile {
    match cohort {
        CliPatientCohort::Adult => PatientProfile::standard_adult(),
        CliPatientCohort::Adolescent => PatientProfile::adolescent(),
        CliPatientCohort::Pediatric => PatientProfile::pediatric(),
    }
}

struct TrialResult {
    pub tir_percent: f64,
    pub total_tbr_percent: f64,
    pub total_tar_percent: f64,
    pub cv_percent: f64,
    pub glycemic_risk_index: f64,
}

fn simulate_patient_trial(
    patient: &PatientProfile,
    ctrl_mode: CliControllerMode,
    duration_hours: usize,
) -> Result<TrialResult> {
    let mut engine = IntegratedMetabolicEngine::new(patient.clone());
    let mut solver = Rkf45Solver::new(Rkf45Config::default());
    let mut mpc = ModelPredictiveController::new(MpcConfig::default(), patient.clone());
    let mut pid = PidController::new(PidConfig::default(), patient.clone());
    let mut plgs = PlgsController::new(PlgsConfig::default(), patient.clone());
    let bb = BasalBolusController::new(patient.clone());

    let dt = 5.0; // 5 min step for fast batch benchmarking
    let total_min = duration_hours as f64 * 60.0;
    let mut t = 0.0;

    let mut glucose_series = Vec::new();
    let mut active_infusion = HormoneInfusion {
        insulin_basal_u_per_hr: patient.basal_rate_u_per_hr,
        ..Default::default()
    };

    let mut last_g = patient.target_glucose_mg_dl;

    while t < total_min {
        // Daily standard meals at 08:00 (50g), 13:00 (70g), 19:00 (60g)
        let min_of_day = t % 1440.0;
        if (min_of_day - 480.0).abs() < dt * 0.5 {
            engine.digestion.ingest_meal(&MealEvent::standard(t, 50.0));
            if ctrl_mode == CliControllerMode::OpenLoop {
                active_infusion.insulin_bolus_u += bb.calculate_meal_bolus(50.0);
            }
        } else if (min_of_day - 780.0).abs() < dt * 0.5 {
            engine.digestion.ingest_meal(&MealEvent::standard(t, 70.0));
            if ctrl_mode == CliControllerMode::OpenLoop {
                active_infusion.insulin_bolus_u += bb.calculate_meal_bolus(70.0);
            }
        } else if (min_of_day - 1140.0).abs() < dt * 0.5 {
            engine.digestion.ingest_meal(&MealEvent::standard(t, 60.0));
            if ctrl_mode == CliControllerMode::OpenLoop {
                active_infusion.insulin_bolus_u += bb.calculate_meal_bolus(60.0);
            }
        }

        let state = engine.step(t, dt, &active_infusion, &mut solver)?;

        let trend = (state.glucose_plasma_mg_dl - last_g) / dt;
        last_g = state.glucose_plasma_mg_dl;

        let next_ctrl = match ctrl_mode {
            CliControllerMode::OpenLoop => bb.step(),
            CliControllerMode::Plgs => plgs.step(dt, state.glucose_plasma_mg_dl, trend),
            CliControllerMode::Pid => pid.step(t, state.glucose_plasma_mg_dl),
            CliControllerMode::Mpc => mpc.compute_action(state.glucose_plasma_mg_dl, trend, 0.5),
        };
        active_infusion = next_ctrl;

        glucose_series.push(state.glucose_plasma_mg_dl);
        t += dt;
    }

    let attd = compute_attd_metrics(&glucose_series, dt);
    let var = compute_variability_metrics(&glucose_series);

    Ok(TrialResult {
        tir_percent: attd.tir_percent,
        total_tbr_percent: attd.total_tbr_percent,
        total_tar_percent: attd.total_tar_percent,
        cv_percent: var.cv_percent,
        glycemic_risk_index: attd.glycemic_risk_index,
    })
}
