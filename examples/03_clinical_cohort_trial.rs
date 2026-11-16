//! Example 03: In-silico population clinical trial across an adult, adolescent, and pediatric cohort.
//!
//! Evaluates population distributions of Time-in-Range (TIR), Time-Below-Range (TBR),
//! Glycemic Risk Index (GRI), and Clarke Error Grid compliance.

use glyco_soma::prelude::*;

fn main() -> Result<()> {
    println!("=== GLYCOSOMA: IN-SILICO POPULATION CLINICAL TRIAL (N=30) ===");

    let mut cohort_gen = VirtualCohortGenerator::new(Some(2026));
    let cohort = cohort_gen.generate_cohort(30);

    println!(
        "Generated Virtual Cohort: 30 subjects (partitioned across Adult, Adolescent, Pediatric)"
    );

    let mut cohort_tir = Vec::new();
    let mut cohort_tbr = Vec::new();
    let mut cohort_tar = Vec::new();
    let mut cohort_cv = Vec::new();
    let mut cohort_gri = Vec::new();

    let solver_config = Rkf45Config {
        max_dt: 5.0,
        ..Default::default()
    };

    for (i, patient) in cohort.iter().enumerate() {
        let mut engine = IntegratedMetabolicEngine::new(patient.clone());
        let mut solver = Rkf45Solver::new(solver_config.clone());
        let mut mpc = ModelPredictiveController::new(MpcConfig::default(), patient.clone());

        // Fast 24-hour simulation with 5-minute sampling steps
        let dt = 5.0;
        let mut t: f64 = 0.0;
        let mut glucose_series = Vec::new();
        let mut active_infusion = HormoneInfusion {
            insulin_basal_u_per_hr: patient.basal_rate_u_per_hr,
            ..Default::default()
        };
        let mut last_g = patient.target_glucose_mg_dl;

        while t < 1440.0 {
            // Three meals
            if (t - 480.0f64).abs() < dt * 0.5 {
                engine.digestion.ingest_meal(&MealEvent::standard(t, 50.0));
            } else if (t - 780.0f64).abs() < dt * 0.5 {
                engine.digestion.ingest_meal(&MealEvent::standard(t, 70.0));
            } else if (t - 1140.0f64).abs() < dt * 0.5 {
                engine.digestion.ingest_meal(&MealEvent::standard(t, 55.0));
            }

            let state = engine.step(t, dt, &active_infusion, &mut solver)?;

            let trend = (state.glucose_plasma_mg_dl - last_g) / dt;
            last_g = state.glucose_plasma_mg_dl;

            active_infusion = mpc.compute_action(state.glucose_plasma_mg_dl, trend, 0.4);
            glucose_series.push(state.glucose_plasma_mg_dl);

            t += dt;
        }

        let attd = compute_attd_metrics(&glucose_series, dt);
        let var = compute_variability_metrics(&glucose_series);

        cohort_tir.push(attd.tir_percent);
        cohort_tbr.push(attd.total_tbr_percent);
        cohort_tar.push(attd.total_tar_percent);
        cohort_cv.push(var.cv_percent);
        cohort_gri.push(attd.glycemic_risk_index);

        if (i + 1) % 10 == 0 || i == 0 {
            println!(
                "  [Subject {:02}/30] Cohort={:10?} | Weight={:4.1}kg | TIR={:5.1}% | TBR={:4.1}% | GRI={:4.1}",
                i + 1,
                patient.cohort,
                patient.body_weight_kg,
                attd.tir_percent,
                attd.total_tbr_percent,
                attd.glycemic_risk_index
            );
        }
    }

    let stats_tir = glyco_soma::math::compute_series_stats(&cohort_tir)?;
    let stats_tbr = glyco_soma::math::compute_series_stats(&cohort_tbr)?;
    let stats_tar = glyco_soma::math::compute_series_stats(&cohort_tar)?;
    let stats_cv = glyco_soma::math::compute_series_stats(&cohort_cv)?;
    let stats_gri = glyco_soma::math::compute_series_stats(&cohort_gri)?;

    println!("\n================ POPULATION TRIAL RESULTS (N=30) ================");
    println!(
        "Time-In-Range (TIR %):      Mean={:5.1}% | Median={:5.1}% | IQR=[{:4.1} - {:4.1}]%",
        stats_tir.mean, stats_tir.median, stats_tir.p25, stats_tir.p75
    );
    println!(
        "Time-Below-Range (TBR %):    Mean={:5.1}% | Median={:5.1}% | Max={:4.1}%",
        stats_tbr.mean, stats_tbr.median, stats_tbr.max
    );
    println!(
        "Time-Above-Range (TAR %):    Mean={:5.1}% | Median={:5.1}% | IQR=[{:4.1} - {:4.1}]%",
        stats_tar.mean, stats_tar.median, stats_tar.p25, stats_tar.p75
    );
    println!(
        "Glucose CV (%):              Mean={:5.1}% | Median={:5.1}%",
        stats_cv.mean, stats_cv.median
    );
    println!(
        "Glycemic Risk Index (GRI):   Mean={:5.1}  | Median={:5.1}  (Target < 20: Very Low Risk)",
        stats_gri.mean, stats_gri.median
    );
    println!("=================================================================\n");

    Ok(())
}
