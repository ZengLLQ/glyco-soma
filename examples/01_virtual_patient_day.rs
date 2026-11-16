//! Example 01: Full 24-hour virtual patient simulation under automated closed-loop MPC control.
//!
//! Simulates a standard adult patient undergoing 3 distinct carbohydrate meals (breakfast 45g, lunch 75g, dinner 55g),
//! a 40-minute moderate aerobic exercise session, and nocturnal dawn phenomenon, managed by an adaptive
//! Model Predictive Controller (MPC) with Continuous Glucose Monitoring (CGM).

use glyco_soma::prelude::*;

fn main() -> Result<()> {
    println!("=== GLYCOSOMA: 24-HOUR CLOSED-LOOP ARTIFICIAL PANCREAS SIMULATION ===");

    let patient = PatientProfile::standard_adult();
    println!("Patient ID:          {}", patient.id);
    println!("Cohort:              {:?}", patient.cohort);
    println!("Body Weight:         {:.1} kg", patient.body_weight_kg);
    println!(
        "Basal Delivery Rate: {:.2} U/hr",
        patient.basal_rate_u_per_hr
    );
    println!("Carb Ratio (CIR):    {:.1} g/U", patient.carb_ratio_g_per_u);
    println!(
        "Sensitivity (ISF):   {:.1} mg/dL/U",
        patient.isf_mg_dl_per_u
    );
    println!(
        "Target Glucose:      {:.1} mg/dL\n",
        patient.target_glucose_mg_dl
    );

    let mut engine = IntegratedMetabolicEngine::new(patient.clone());
    let mut cgm = CgmSensor::new(
        CgmParams::default(),
        patient.target_glucose_mg_dl,
        Some(101),
    );
    let mut solver = Rkf45Solver::new(Rkf45Config::default());
    let mut mpc = ModelPredictiveController::new(MpcConfig::default(), patient.clone());

    // Schedule meals: Breakfast at 07:30 (45g), Lunch at 12:30 (75g), Dinner at 19:00 (55g)
    let meals = [
        MealEvent::standard(7.5 * 60.0, 45.0),
        MealEvent::slow_high_fat(12.5 * 60.0, 75.0),
        MealEvent::standard(19.0 * 60.0, 55.0),
    ];

    // Schedule 40-min aerobic exercise at 16:00
    engine
        .exercise
        .trigger_exercise(ExerciseEvent::moderate_aerobic(16.0 * 60.0, 40.0));

    let total_time_min = 24.0 * 60.0;
    let dt = 1.0; // 1-minute step

    let mut time_series = Vec::new();
    let mut plasma_glucose_series = Vec::new();
    let mut cgm_series = Vec::new();
    let mut ref_series = Vec::new();

    let mut t = 0.0;
    let mut active_infusion = HormoneInfusion {
        insulin_basal_u_per_hr: patient.basal_rate_u_per_hr,
        ..Default::default()
    };

    let mut last_cgm = patient.target_glucose_mg_dl;
    let mut last_cgm_time = 0.0;

    while t < total_time_min {
        // Ingest meals
        for meal in &meals {
            if (t - meal.start_time_min).abs() < dt * 0.5 {
                engine.digestion.ingest_meal(meal);
                println!(
                    "[{:02}:{:02}] Meal ingested: {:.0}g carbs (GI: {:.2})",
                    (t / 60.0).floor() as u32,
                    (t % 60.0).floor() as u32,
                    meal.carbs_grams,
                    meal.glycemic_index
                );
            }
        }

        // Exercise announcement
        if (t - 16.0 * 60.0).abs() < dt * 0.5 {
            println!(
                "[{:02}:{:02}] Exercise started: 40 min aerobic cycling (5 METs)",
                (t / 60.0).floor() as u32,
                (t % 60.0).floor() as u32
            );
        }

        // Advance physiology
        let state = engine.step(t, dt, &active_infusion, &mut solver)?;

        // Clear discrete bolus commands
        active_infusion.insulin_bolus_u = 0.0;
        active_infusion.glucagon_bolus_ug = 0.0;

        // Sample CGM sensor
        if let Some(cgm_reading) = cgm.step(t, dt, state.glucose_plasma_mg_dl) {
            let trend = (cgm_reading - last_cgm) / (t - last_cgm_time).max(1.0);
            last_cgm = cgm_reading;
            last_cgm_time = t;

            // MPC controller calculates optimal next infusion
            active_infusion = mpc.compute_action(cgm_reading, trend, 0.4);

            cgm_series.push(cgm_reading);
            ref_series.push(state.glucose_plasma_mg_dl);
        }

        time_series.push(t);
        plasma_glucose_series.push(state.glucose_plasma_mg_dl);

        t += dt;
    }

    // Render ASCII strip chart
    let chart =
        render_ascii_glucose_chart(&time_series, &plasma_glucose_series, 74, 18, 70.0, 180.0);
    println!("\n{}", chart);

    // Compute Clinical Metrics
    let attd = compute_attd_metrics(&plasma_glucose_series, dt);
    let var = compute_variability_metrics(&plasma_glucose_series);
    let clark = evaluate_clark_grid(&ref_series, &cgm_series);

    println!("\n--- 24-HOUR CLINICAL OUTCOME REPORT ---");
    println!(
        "Time In Range (70-180 mg/dL):  {:5.1}% (Target > 70%)",
        attd.tir_percent
    );
    println!(
        "Time Below Range (< 70 mg/dL): {:5.1}% (Target < 4%)",
        attd.total_tbr_percent
    );
    println!(
        "Time Above Range (> 180 mg/dL):{:5.1}% (Target < 25%)",
        attd.total_tar_percent
    );
    println!(
        "Mean Glucose:                  {:5.1} mg/dL",
        attd.mean_glucose_mg_dl
    );
    println!("Estimated HbA1c (GMI):         {:4.1}%", attd.gmi_percent);
    println!(
        "Glycemic Variability (CV%):    {:5.1}% (Target <= 36%)",
        var.cv_percent
    );
    println!(
        "Glycemic Risk Index (GRI):     {:4.1}/100 [{}]",
        attd.glycemic_risk_index, attd.gri_category
    );
    println!(
        "Clarke Grid Clinical Accuracy: {:5.1}% in Zone A+B (Target >= 95%)",
        clark.percent_a_plus_b
    );
    println!("---------------------------------------\n");

    Ok(())
}
