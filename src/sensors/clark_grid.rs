//! Clarke Error Grid Analysis (CEGA) and Consensus Error Grid clinical accuracy validation.
//!
//! # Zones
//! - **Zone A**: Values within 20% of reference, or both $\le 70$ mg/dL (clinically accurate, no treatment error).
//! - **Zone B**: Values outside 20% error but would lead to benign or no treatment error.
//! - **Zone C**: Values that would result in unnecessary or overcorrecting clinical action.
//! - **Zone D**: Values that fail to detect clinically significant hypoglycemia or hyperglycemia.
//! - **Zone E**: Values that confuse hypoglycemia for hyperglycemia (or vice versa), leading to dangerous treatment inversion.

use serde::{Deserialize, Serialize};

/// Clarke Error Grid zone classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ClarkZone {
    ZoneA,
    ZoneB,
    ZoneC,
    ZoneD,
    ZoneE,
}

impl ClarkZone {
    /// Text label.
    pub fn label(self) -> &'static str {
        match self {
            ClarkZone::ZoneA => "Zone A (Accurate)",
            ClarkZone::ZoneB => "Zone B (Benign)",
            ClarkZone::ZoneC => "Zone C (Overcorrection)",
            ClarkZone::ZoneD => "Zone D (Dangerous Failure to Detect)",
            ClarkZone::ZoneE => "Zone E (Erroneous Treatment)",
        }
    }

    /// True if clinically acceptable (Zone A or B).
    pub fn is_clinically_acceptable(self) -> bool {
        matches!(self, ClarkZone::ZoneA | ClarkZone::ZoneB)
    }
}

/// Aggregated Clarke Error Grid evaluation summary.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ClarkGridReport {
    /// Total data pairs evaluated.
    pub total_pairs: usize,
    /// Count in Zone A.
    pub count_a: usize,
    /// Count in Zone B.
    pub count_b: usize,
    /// Count in Zone C.
    pub count_c: usize,
    /// Count in Zone D.
    pub count_d: usize,
    /// Count in Zone E.
    pub count_e: usize,
    /// Percentage in Zone A.
    pub percent_a: f64,
    /// Percentage in Zone B.
    pub percent_b: f64,
    /// Percentage in Zone A + B (ISO 15197 target: $\ge 95\%$).
    pub percent_a_plus_b: f64,
    /// Percentage in Zone C.
    pub percent_c: f64,
    /// Percentage in Zone D.
    pub percent_d: f64,
    /// Percentage in Zone E.
    pub percent_e: f64,
    /// ISO 15197 clinical compliance verdict (Zone A+B >= 95% and Zone E == 0).
    pub iso_compliant: bool,
}

/// Classify a single pair $(G_{ref}, G_{meas})$ into its Clarke Error Grid zone.
///
/// Both inputs must be in mg/dL.
pub fn classify_clark_point(ref_val: f64, meas_val: f64) -> ClarkZone {
    let r = ref_val.max(1.0);
    let m = meas_val.max(1.0);

    // 1. Zone A boundary:
    // (a) Within 20% of reference: |m - r| <= 0.20 * r
    // (b) Or both reference and measurement <= 70 mg/dL: |m - r| <= 14 mg/dL
    if (m >= 0.80 * r && m <= 1.20 * r) || (r <= 70.0 && m <= 70.0) {
        return ClarkZone::ZoneA;
    }

    // 2. Zone E boundaries (most severe inversion):
    // Reference <= 70 and measurement >= 180 (hypo treated as hyper)
    if r <= 70.0 && m >= 180.0 {
        return ClarkZone::ZoneE;
    }
    // Reference >= 180 and measurement <= 70 (hyper treated as hypo)
    if r >= 180.0 && m <= 70.0 {
        return ClarkZone::ZoneE;
    }

    // 3. Zone D boundaries (dangerous failure to detect):
    // Reference <= 70 and measurement in (70, 180) -> missed hypoglycemia
    if r <= 70.0 && m > 70.0 && m < 180.0 {
        return ClarkZone::ZoneD;
    }
    // Reference >= 240 and measurement in [70, 180] -> missed severe hyperglycemia
    if r >= 240.0 && (70.0..=180.0).contains(&m) {
        return ClarkZone::ZoneD;
    }

    // 4. Zone C boundaries (overcorrection / unnecessary action):
    // Upper Zone C: Reference between 70 and 290, and measurement > reference + 110
    if (70.0..=290.0).contains(&r) && m > (r + 110.0) {
        return ClarkZone::ZoneC;
    }
    // Lower Zone C: Reference between 130 and 180, and measurement < 70 / 0.8
    if (130.0..=180.0).contains(&r) && m < (7.0 / 5.0) * (r - 130.0) {
        return ClarkZone::ZoneC;
    }

    // 5. Zone B: Remaining points fall into the benign Zone B
    ClarkZone::ZoneB
}

/// Evaluate an entire time series of reference and measured glucose pairs.
pub fn evaluate_clark_grid(ref_vals: &[f64], meas_vals: &[f64]) -> ClarkGridReport {
    let n = ref_vals.len().min(meas_vals.len());
    if n == 0 {
        return ClarkGridReport::default();
    }

    let mut count_a = 0;
    let mut count_b = 0;
    let mut count_c = 0;
    let mut count_d = 0;
    let mut count_e = 0;

    for i in 0..n {
        match classify_clark_point(ref_vals[i], meas_vals[i]) {
            ClarkZone::ZoneA => count_a += 1,
            ClarkZone::ZoneB => count_b += 1,
            ClarkZone::ZoneC => count_c += 1,
            ClarkZone::ZoneD => count_d += 1,
            ClarkZone::ZoneE => count_e += 1,
        }
    }

    let total = n as f64;
    let percent_a = (count_a as f64 / total) * 100.0;
    let percent_b = (count_b as f64 / total) * 100.0;
    let percent_c = (count_c as f64 / total) * 100.0;
    let percent_d = (count_d as f64 / total) * 100.0;
    let percent_e = (count_e as f64 / total) * 100.0;
    let percent_a_plus_b = percent_a + percent_b;

    let iso_compliant = percent_a_plus_b >= 95.0 && count_e == 0;

    ClarkGridReport {
        total_pairs: n,
        count_a,
        count_b,
        count_c,
        count_d,
        count_e,
        percent_a,
        percent_b,
        percent_a_plus_b,
        percent_c,
        percent_d,
        percent_e,
        iso_compliant,
    }
}
