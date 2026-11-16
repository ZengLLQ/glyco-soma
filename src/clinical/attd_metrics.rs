//! Advanced Technologies & Treatments for Diabetes (ATTD 2024 / ADA) consensus clinical metrics.
//!
//! # Standard Clinical Targets
//! - **Time In Range (TIR)**: $70 - 180$ mg/dL ($> 70\%$ target).
//! - **Time Below Range Level 1 (TBR L1)**: $54 - 69$ mg/dL ($< 4\%$ target).
//! - **Time Below Range Level 2 (TBR L2)**: $< 54$ mg/dL ($< 1\%$ target).
//! - **Time Above Range Level 1 (TAR L1)**: $181 - 250$ mg/dL ($< 25\%$ target).
//! - **Time Above Range Level 2 (TAR L2)**: $> 250$ mg/dL ($< 5\%$ target).
//! - **Glucose Management Indicator (GMI)**: $GMI = 3.31 + 0.02392 \times \bar{G}$ (estimated HbA1c).
//! - **Glycemic Risk Index (GRI)**: $3.0 \times TBR_{L2} + 2.4 \times TBR_{L1} + 1.6 \times TAR_{L2} + 0.8 \times TAR_{L1}$.

use serde::{Deserialize, Serialize};

/// Comprehensive ATTD / ADA consensus report.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct AttdMetricsReport {
    /// Total duration of recorded sensor data in hours.
    pub total_hours: f64,
    /// Total number of sensor readings analyzed.
    pub sample_count: usize,
    /// Mean glucose concentration in mg/dL.
    pub mean_glucose_mg_dl: f64,
    /// Median glucose concentration in mg/dL.
    pub median_glucose_mg_dl: f64,
    /// Glucose Management Indicator (GMI, estimated % HbA1c).
    pub gmi_percent: f64,
    /// Time in Range (TIR, 70-180 mg/dL) percentage.
    pub tir_percent: f64,
    /// Time Below Range Level 1 (TBR L1, 54-69 mg/dL) percentage.
    pub tbr_l1_percent: f64,
    /// Time Below Range Level 2 (TBR L2, < 54 mg/dL) percentage.
    pub tbr_l2_percent: f64,
    /// Total Time Below Range (TBR = L1 + L2) percentage.
    pub total_tbr_percent: f64,
    /// Time Above Range Level 1 (TAR L1, 181-250 mg/dL) percentage.
    pub tar_l1_percent: f64,
    /// Time Above Range Level 2 (TAR L2, > 250 mg/dL) percentage.
    pub tar_l2_percent: f64,
    /// Total Time Above Range (TAR = L1 + L2) percentage.
    pub total_tar_percent: f64,
    /// Glycemic Risk Index (GRI) score (0 to 100).
    pub glycemic_risk_index: f64,
    /// Hypoglycemia Component of GRI (CHypo = 3.0 * TBR_L2 + 2.4 * TBR_L1).
    pub c_hypo: f64,
    /// Hyperglycemia Component of GRI (CHyper = 1.6 * TAR_L2 + 0.8 * TAR_L1).
    pub c_hyper: f64,
    /// GRI clinical risk percentile category.
    pub gri_category: String,
    /// Clinical consensus compliance (TIR >= 70% and Total TBR < 4%).
    pub consensus_compliant: bool,
}

/// Compute ATTD consensus metrics from a slice of glucose readings (in mg/dL) sampled at fixed interval `sampling_min`.
pub fn compute_attd_metrics(glucose_series: &[f64], sampling_min: f64) -> AttdMetricsReport {
    if glucose_series.is_empty() {
        return AttdMetricsReport::default();
    }

    let n = glucose_series.len();
    let sum: f64 = glucose_series.iter().sum();
    let mean_glucose = sum / n as f64;

    let mut sorted = glucose_series.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_glucose = sorted[n / 2];

    // GMI formula (Bergenstal et al., Diabetes Care 2018): GMI(%) = 3.31 + 0.02392 * mean_glucose
    let gmi_percent = 3.31 + 0.02392 * mean_glucose;

    let mut count_tbr_l2 = 0; // < 54
    let mut count_tbr_l1 = 0; // 54 - 69
    let mut count_tir = 0; // 70 - 180
    let mut count_tar_l1 = 0; // 181 - 250
    let mut count_tar_l2 = 0; // > 250

    for &g in glucose_series {
        if g < 54.0 {
            count_tbr_l2 += 1;
        } else if g < 70.0 {
            count_tbr_l1 += 1;
        } else if g <= 180.0 {
            count_tir += 1;
        } else if g <= 250.0 {
            count_tar_l1 += 1;
        } else {
            count_tar_l2 += 1;
        }
    }

    let total = n as f64;
    let tbr_l2_percent = (count_tbr_l2 as f64 / total) * 100.0;
    let tbr_l1_percent = (count_tbr_l1 as f64 / total) * 100.0;
    let total_tbr_percent = tbr_l1_percent + tbr_l2_percent;

    let tir_percent = (count_tir as f64 / total) * 100.0;

    let tar_l1_percent = (count_tar_l1 as f64 / total) * 100.0;
    let tar_l2_percent = (count_tar_l2 as f64 / total) * 100.0;
    let total_tar_percent = tar_l1_percent + tar_l2_percent;

    // Glycemic Risk Index (Klonoff et al., J Diabetes Sci Technol 2022)
    // GRI = (3.0 * TBR_L2) + (2.4 * TBR_L1) + (1.6 * TAR_L2) + (0.8 * TAR_L1)
    let c_hypo = (3.0 * tbr_l2_percent) + (2.4 * tbr_l1_percent);
    let c_hyper = (1.6 * tar_l2_percent) + (0.8 * tar_l1_percent);
    let raw_gri = c_hypo + c_hyper;
    let glycemic_risk_index = raw_gri.clamp(0.0, 100.0);

    let gri_category = if glycemic_risk_index < 20.0 {
        "Very Low Risk (Zone A)".to_string()
    } else if glycemic_risk_index < 40.0 {
        "Low Risk (Zone B)".to_string()
    } else if glycemic_risk_index < 60.0 {
        "Moderate Risk (Zone C)".to_string()
    } else if glycemic_risk_index < 80.0 {
        "High Risk (Zone D)".to_string()
    } else {
        "Very High Risk (Zone E)".to_string()
    };

    let total_hours = (n as f64 * sampling_min) / 60.0;
    let consensus_compliant =
        tir_percent >= 70.0 && total_tbr_percent <= 4.0 && tbr_l2_percent <= 1.0;

    AttdMetricsReport {
        total_hours,
        sample_count: n,
        mean_glucose_mg_dl: mean_glucose,
        median_glucose_mg_dl: median_glucose,
        gmi_percent,
        tir_percent,
        tbr_l1_percent,
        tbr_l2_percent,
        total_tbr_percent,
        tar_l1_percent,
        tar_l2_percent,
        total_tar_percent,
        glycemic_risk_index,
        c_hypo,
        c_hyper,
        gri_category,
        consensus_compliant,
    }
}
