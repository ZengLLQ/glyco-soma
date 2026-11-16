//! Multi-day full closed-loop simulation integration tests.

use glyco_soma::prelude::*;

#[test]
fn test_multiday_closed_loop_mpc_stability() {
    let patient = PatientProfile::standard_adult();
    let mut engine = IntegratedMetabolicEngine::new(patient.clone());
    let mut cgm = CgmSensor::new(CgmParams::default(), patient.target_glucose_mg_dl, Some(77));
    let mut solver = Rkf45Solver::new(Rkf45Config::default());
    let mut mpc = ModelPredictiveController::new(MpcConfig::default(), patient.clone());

    // 48-hour simulation (2 days)
    let total_time_min = 48.0 * 60.0;
    let dt = 5.0; // 5 min steps
    let mut t = 0.0;

    let mut glucose_history = Vec::new();
    let mut active_infusion = HormoneInfusion {
        insulin_basal_u_per_hr: patient.basal_rate_u_per_hr,
        ..Default::default()
    };

    let mut last_cgm = patient.target_glucose_mg_dl;
    let mut last_time = 0.0;
    let mut iob = 0.0;

    while t < total_time_min {
        // Daily meals at 08:00 (50g), 13:00 (70g), 19:00 (60g) with prandial bolus (55% closed-loop split)
        let min_of_day: f64 = t % 1440.0;
        let mut meal_bolus = 0.0;
        if (min_of_day - 480.0f64).abs() < dt * 0.5 {
            engine.digestion.ingest_meal(&MealEvent::standard(t, 50.0));
            meal_bolus = (50.0 / patient.carb_ratio_g_per_u) * 0.55;
        } else if (min_of_day - 780.0f64).abs() < dt * 0.5 {
            engine.digestion.ingest_meal(&MealEvent::standard(t, 70.0));
            meal_bolus = (70.0 / patient.carb_ratio_g_per_u) * 0.55;
        } else if (min_of_day - 1140.0f64).abs() < dt * 0.5 {
            engine.digestion.ingest_meal(&MealEvent::standard(t, 60.0));
            meal_bolus = (60.0 / patient.carb_ratio_g_per_u) * 0.55;
        }
        active_infusion.insulin_bolus_u = meal_bolus;
        iob = iob * (-dt / 65.0).exp() + meal_bolus;

        let state = engine.step(t, dt, &active_infusion, &mut solver).unwrap();

        if let Some(reading) = cgm.step(t, dt, state.glucose_plasma_mg_dl) {
            let trend = (reading - last_cgm) / (t - last_time).max(1.0);
            last_cgm = reading;
            last_time = t;

            active_infusion = mpc.compute_action(reading, trend, iob);
        }

        glucose_history.push(state.glucose_plasma_mg_dl);
        t += dt;
    }

    let attd = compute_attd_metrics(&glucose_history, dt);
    let var = compute_variability_metrics(&glucose_history);

    // Assert clinical targets
    assert!(
        attd.tir_percent > 70.0,
        "TIR {:.1}% must exceed 70% consensus target",
        attd.tir_percent
    );
    assert!(
        attd.total_tbr_percent < 5.0,
        "Total TBR {:.1}% must remain below 5%",
        attd.total_tbr_percent
    );
    assert!(
        var.cv_percent < 36.0,
        "Glucose CV {:.1}% must satisfy stability <= 36%",
        var.cv_percent
    );
}

#[test]
fn test_exercise_glycemic_drop_and_recovery() {
    let patient = PatientProfile::standard_adult();
    let mut engine = IntegratedMetabolicEngine::new(patient.clone());
    let mut solver = Rkf45Solver::new(Rkf45Config::default());

    // Schedule 60-min intense exercise (8 METs) at t = 60 min
    engine
        .exercise
        .trigger_exercise(ExerciseEvent::vigorous(60.0, 60.0));

    let infusion = HormoneInfusion {
        insulin_basal_u_per_hr: patient.basal_rate_u_per_hr,
        ..Default::default()
    };

    let dt = 1.0;
    let mut t = 0.0;
    let mut min_exercise_glucose = f64::MAX;

    while t < 240.0 {
        let state = engine.step(t, dt, &infusion, &mut solver).unwrap();
        if (60.0..=120.0).contains(&t) && state.glucose_plasma_mg_dl < min_exercise_glucose {
            min_exercise_glucose = state.glucose_plasma_mg_dl;
        }
        t += dt;
    }

    // Vigorous exercise without carb intake should drop glucose
    assert!(
        min_exercise_glucose < patient.target_glucose_mg_dl,
        "Exercise must induce non-insulin mediated glucose uptake"
    );
}
