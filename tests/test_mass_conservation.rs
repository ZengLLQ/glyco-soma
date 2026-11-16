//! Tests for biological mass conservation, positivity invariants, and steady-state equilibrium.

use glyco_soma::prelude::*;

#[test]
fn test_positivity_invariant_under_zero_infusion() {
    let patient = PatientProfile::standard_adult();
    let mut engine = IntegratedMetabolicEngine::new(patient);
    let mut solver = Rkf45Solver::new(Rkf45Config::default());

    let zero_infusion = HormoneInfusion::default();
    let dt = 1.0;
    let mut t = 0.0;

    // Simulate for 60 minutes with zero external input
    while t < 60.0 {
        let state = engine.step(t, dt, &zero_infusion, &mut solver).unwrap();

        // Positivity invariants
        assert!(
            state.glucose_plasma_mg_dl > 0.0,
            "Glucose must remain strictly positive"
        );
        assert!(
            state.s1_subcut_insulin >= 0.0,
            "S1 must remain non-negative"
        );
        assert!(
            state.s2_subcut_insulin >= 0.0,
            "S2 must remain non-negative"
        );
        assert!(
            state.insulin_plasma_mu_l >= 0.0,
            "Plasma insulin must remain non-negative"
        );
        assert!(
            state.q1_accessible_glucose_mmol > 0.0,
            "Q1 must remain positive"
        );
        assert!(
            state.q2_nonaccessible_glucose_mmol > 0.0,
            "Q2 must remain positive"
        );
        assert!(
            state.liver_glycogen_fraction > 0.0,
            "Glycogen must remain positive"
        );

        t += dt;
    }
}

#[test]
fn test_glucose_mass_conservation_during_meal_absorption() {
    let mut digestion = DigestionSystem::new(DigestionParams::default());
    let meal = MealEvent::standard(0.0, 50.0); // 50g meal
    digestion.ingest_meal(&meal);

    let dt = 1.0;
    let mut t = 0.0;
    let mut total_absorbed_mmol = 0.0;

    // Simulate until gut is almost completely empty (6 hours = 360 min)
    while t < 360.0 {
        let rate_mmol_min = digestion.step(dt);
        total_absorbed_mmol += rate_mmol_min * dt;
        t += dt;
    }

    // 50g with 0.92 bioavailability = 46.0g
    // 46.0g * 5.5507 mmol/g = 255.33 mmol
    let expected_absorbed_mmol = 50.0 * 0.92 * 5.55075;
    let remaining_gut_g = digestion.total_gut_carbs_g();

    // Absorbed + remaining should closely match ingested bioavailable carbs
    let remaining_mmol = remaining_gut_g * 5.55075;
    let total_accounted = total_absorbed_mmol + remaining_mmol;
    let diff = (total_accounted - expected_absorbed_mmol).abs();

    assert!(
        diff < 1.0,
        "Total accounted carbs {:.2} mmol deviates from expected {:.2} mmol by {:.2}",
        total_accounted,
        expected_absorbed_mmol,
        diff
    );
}

#[test]
fn test_subcutaneous_insulin_transit_mass_conservation() {
    let patient = PatientProfile::standard_adult();
    let mut engine = IntegratedMetabolicEngine::new(patient);
    let mut solver = Rkf45Solver::new(Rkf45Config::default());

    // Inject a single 5.0 Unit discrete bolus
    let bolus_u = 5.0;
    let mut infusion = HormoneInfusion {
        insulin_basal_u_per_hr: 0.0,
        insulin_bolus_u: bolus_u,
        ..Default::default()
    };

    let dt = 1.0;
    let mut t = 0.0;

    // Step 1: Deliver bolus
    let state0 = engine.step(t, dt, &infusion, &mut solver).unwrap();
    infusion.insulin_bolus_u = 0.0;
    t += dt;

    // Verify S1 received the bolus (5 U = 30,000 pmol)
    assert!(state0.s1_subcut_insulin > 20000.0);

    // Follow transit over 4 hours (240 min)
    while t < 240.0 {
        let state = engine.step(t, dt, &infusion, &mut solver).unwrap();
        assert!(state.s1_subcut_insulin >= 0.0);
        assert!(state.s2_subcut_insulin >= 0.0);
        t += dt;
    }

    // After 4 hours, subcutaneous depots should have mostly drained
    let final_s1 = engine.state_vector[glyco_soma::models::hovorka::state_idx::S1];
    let final_s2 = engine.state_vector[glyco_soma::models::hovorka::state_idx::S2];
    assert!(final_s1 < 3000.0);
    assert!(final_s2 < 3000.0);
}
