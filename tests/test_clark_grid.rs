//! Verification tests for Clarke Error Grid zone boundary classifications.

use glyco_soma::sensors::clark_grid::{classify_clark_point, evaluate_clark_grid, ClarkZone};

#[test]
fn test_zone_a_exact_and_within_20_percent() {
    // Exact match
    assert_eq!(classify_clark_point(100.0, 100.0), ClarkZone::ZoneA);
    // +15% error
    assert_eq!(classify_clark_point(100.0, 115.0), ClarkZone::ZoneA);
    // -15% error
    assert_eq!(classify_clark_point(100.0, 85.0), ClarkZone::ZoneA);
    // Boundary 20%
    assert_eq!(classify_clark_point(100.0, 120.0), ClarkZone::ZoneA);
    assert_eq!(classify_clark_point(100.0, 80.0), ClarkZone::ZoneA);

    // Both <= 70 mg/dL (hypoglycemia tolerance)
    assert_eq!(classify_clark_point(50.0, 65.0), ClarkZone::ZoneA);
    assert_eq!(classify_clark_point(65.0, 50.0), ClarkZone::ZoneA);
}

#[test]
fn test_zone_e_dangerous_treatment_inversion() {
    // Reference in severe hypo, sensor reports severe hyper (insulin given to crashing patient)
    assert_eq!(classify_clark_point(55.0, 220.0), ClarkZone::ZoneE);
    assert_eq!(classify_clark_point(68.0, 185.0), ClarkZone::ZoneE);

    // Reference in hyper, sensor reports hypo (sugar given to hyperglycemic patient)
    assert_eq!(classify_clark_point(250.0, 50.0), ClarkZone::ZoneE);
    assert_eq!(classify_clark_point(190.0, 65.0), ClarkZone::ZoneE);
}

#[test]
fn test_zone_d_failure_to_detect() {
    // Missed hypoglycemia (ref is hypo <= 70, sensor reads normal 70-180)
    assert_eq!(classify_clark_point(50.0, 110.0), ClarkZone::ZoneD);
    assert_eq!(classify_clark_point(65.0, 140.0), ClarkZone::ZoneD);

    // Missed severe hyperglycemia (ref is hyper >= 240, sensor reads normal 70-180)
    assert_eq!(classify_clark_point(280.0, 120.0), ClarkZone::ZoneD);
    assert_eq!(classify_clark_point(320.0, 150.0), ClarkZone::ZoneD);
}

#[test]
fn test_zone_c_overcorrection() {
    // Overcorrecting treatment (sensor reads much higher than reference)
    assert_eq!(classify_clark_point(120.0, 250.0), ClarkZone::ZoneC);
    assert_eq!(classify_clark_point(150.0, 280.0), ClarkZone::ZoneC);
}

#[test]
fn test_zone_b_benign_errors() {
    // Outside 20% but within benign decision limits
    assert_eq!(classify_clark_point(100.0, 130.0), ClarkZone::ZoneB);
    assert_eq!(classify_clark_point(150.0, 110.0), ClarkZone::ZoneB);
    assert_eq!(classify_clark_point(200.0, 250.0), ClarkZone::ZoneB);
}

#[test]
fn test_clark_grid_aggregate_report() {
    let ref_vals = vec![100.0, 110.0, 55.0, 220.0, 130.0, 100.0];
    let meas_vals = vec![105.0, 108.0, 58.0, 215.0, 135.0, 102.0];

    let report = evaluate_clark_grid(&ref_vals, &meas_vals);
    assert_eq!(report.total_pairs, 6);
    assert_eq!(report.count_a, 6);
    assert_eq!(report.percent_a, 100.0);
    assert!(report.iso_compliant);
}
