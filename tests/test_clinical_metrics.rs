//! Unit tests for clinical ATTD consensus metrics, GMI, GRI, MAGE, and Kovatchev indices.

use glyco_soma::clinical::attd_metrics::compute_attd_metrics;
use glyco_soma::clinical::variability::compute_variability_metrics;

#[test]
fn test_attd_metrics_perfect_range() {
    // 100 samples all within 100-110 mg/dL (100% in range)
    let series = vec![105.0; 100];
    let report = compute_attd_metrics(&series, 5.0);

    assert_eq!(report.sample_count, 100);
    assert!((report.mean_glucose_mg_dl - 105.0).abs() < 1e-6);
    assert_eq!(report.tir_percent, 100.0);
    assert_eq!(report.total_tbr_percent, 0.0);
    assert_eq!(report.total_tar_percent, 0.0);
    assert!(report.consensus_compliant);

    // GMI = 3.31 + 0.02392 * 105 = 5.82%
    assert!((report.gmi_percent - 5.82).abs() < 0.05);
    assert_eq!(report.glycemic_risk_index, 0.0);
}

#[test]
fn test_attd_metrics_hypo_and_hyper_breakdown() {
    let mut series = Vec::new();
    // 10 severe hypo (< 54)
    series.extend(vec![48.0; 10]);
    // 10 mild hypo (54-69)
    series.extend(vec![62.0; 10]);
    // 60 in range (70-180)
    series.extend(vec![110.0; 60]);
    // 15 elevated (181-250)
    series.extend(vec![200.0; 15]);
    // 5 very high (> 250)
    series.extend(vec![290.0; 5]);

    let report = compute_attd_metrics(&series, 5.0);

    assert_eq!(report.sample_count, 100);
    assert_eq!(report.tbr_l2_percent, 10.0);
    assert_eq!(report.tbr_l1_percent, 10.0);
    assert_eq!(report.total_tbr_percent, 20.0);
    assert_eq!(report.tir_percent, 60.0);
    assert_eq!(report.tar_l1_percent, 15.0);
    assert_eq!(report.tar_l2_percent, 5.0);
    assert_eq!(report.total_tar_percent, 20.0);

    // CHypo = 3.0 * 10 + 2.4 * 10 = 54.0
    // CHyper = 1.6 * 5 + 0.8 * 15 = 8 + 12 = 20.0
    // GRI = 74.0
    assert!((report.c_hypo - 54.0).abs() < 1e-4);
    assert!((report.c_hyper - 20.0).abs() < 1e-4);
    assert!((report.glycemic_risk_index - 74.0).abs() < 1e-4);
    assert!(!report.consensus_compliant);
}

#[test]
fn test_variability_metrics_cv_and_kovatchev() {
    let series = vec![80.0, 95.0, 110.0, 125.0, 140.0, 125.0, 110.0, 95.0];
    let var = compute_variability_metrics(&series);

    assert!(var.cv_percent > 0.0);
    assert!(var.std_dev_mg_dl > 0.0);
    assert!(var.is_stable_cv); // CV < 36%
    assert!(var.lbgi >= 0.0);
    assert!(var.hbgi >= 0.0);
}

#[test]
fn test_kovatchev_high_hypo_risk() {
    // Series with extreme hypoglycemia
    let series = vec![42.0, 45.0, 48.0, 50.0, 52.0];
    let var = compute_variability_metrics(&series);

    // LBGI should be very high (> 5.0)
    assert!(var.lbgi > 5.0);
    assert_eq!(var.lbgi_risk_category, "High Risk (> 5.0)");
}
