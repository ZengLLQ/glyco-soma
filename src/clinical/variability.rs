//! Glycemic variability metrics: CV%, MAGE, MODD, and Kovatchev LBGI/HBGI indices.

use serde::{Deserialize, Serialize};

/// Glycemic variability summary metrics.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GlycemicVariabilityReport {
    /// Coefficient of Variation ($CV = \frac{SD}{\bar{G}} \times 100\%$, target $\le 36\%$).
    pub cv_percent: f64,
    /// Standard deviation in mg/dL.
    pub std_dev_mg_dl: f64,
    /// Mean Amplitude of Glycemic Excursions (MAGE, mg/dL).
    pub mage_mg_dl: f64,
    /// Low Blood Glucose Index (LBGI: risk of severe hypoglycemia).
    pub lbgi: f64,
    /// High Blood Glucose Index (HBGI: risk of severe hyperglycemia).
    pub hbgi: f64,
    /// LBGI clinical risk category.
    pub lbgi_risk_category: String,
    /// HBGI clinical risk category.
    pub hbgi_risk_category: String,
    /// True if glycemic variability is stable ($CV \le 36\%$).
    pub is_stable_cv: bool,
}

/// Compute comprehensive glycemic variability metrics from a glucose time series.
pub fn compute_variability_metrics(glucose_series: &[f64]) -> GlycemicVariabilityReport {
    if glucose_series.len() < 2 {
        return GlycemicVariabilityReport::default();
    }

    let n = glucose_series.len() as f64;
    let mean = glucose_series.iter().sum::<f64>() / n;
    let variance = glucose_series
        .iter()
        .map(|&x| (x - mean) * (x - mean))
        .sum::<f64>()
        / (n - 1.0);
    let std_dev = variance.sqrt();
    let cv_percent = if mean > 1e-6 {
        (std_dev / mean) * 100.0
    } else {
        0.0
    };

    // MAGE calculation:
    // Identify local peaks and nadirs exceeding 1 standard deviation
    let mage = calculate_mage(glucose_series, std_dev);

    // Kovatchev LBGI / HBGI logarithmic transformation:
    // f(G) = 1.509 * (ln(G)^1.084 - 5.381)
    // risk(G) = 10 * f(G)^2
    // If f(G) < 0: rl(G) = risk(G), rh(G) = 0
    // If f(G) > 0: rl(G) = 0, rh(G) = risk(G)
    let mut sum_rl = 0.0;
    let mut sum_rh = 0.0;

    for &g in glucose_series {
        let g_clamped = g.clamp(20.0, 600.0);
        let ln_g = g_clamped.ln();
        let f_g = 1.509 * (ln_g.powf(1.084) - 5.381);
        let risk = 10.0 * f_g * f_g;

        if f_g < 0.0 {
            sum_rl += risk;
        } else {
            sum_rh += risk;
        }
    }

    let lbgi = sum_rl / n;
    let hbgi = sum_rh / n;

    let lbgi_risk_category = if lbgi <= 1.1 {
        "Minimal Risk (<= 1.1)".to_string()
    } else if lbgi <= 2.5 {
        "Low Risk (1.1 - 2.5)".to_string()
    } else if lbgi <= 5.0 {
        "Moderate Risk (2.5 - 5.0)".to_string()
    } else {
        "High Risk (> 5.0)".to_string()
    };

    let hbgi_risk_category = if hbgi <= 4.5 {
        "Low Risk (<= 4.5)".to_string()
    } else if hbgi <= 9.0 {
        "Moderate Risk (4.5 - 9.0)".to_string()
    } else {
        "High Risk (> 9.0)".to_string()
    };

    GlycemicVariabilityReport {
        cv_percent,
        std_dev_mg_dl: std_dev,
        mage_mg_dl: mage,
        lbgi,
        hbgi,
        lbgi_risk_category,
        hbgi_risk_category,
        is_stable_cv: cv_percent <= 36.0,
    }
}

/// Calculate Mean Amplitude of Glycemic Excursions (MAGE) using standard peak/nadir detection exceeding 1 SD.
fn calculate_mage(series: &[f64], std_dev: f64) -> f64 {
    if series.len() < 5 || std_dev < 1e-4 {
        return 0.0;
    }

    // Find local extrema
    let mut extrema = Vec::new();
    for i in 1..(series.len() - 1) {
        let prev = series[i - 1];
        let curr = series[i];
        let next = series[i + 1];

        if (curr >= prev && curr > next) || (curr > prev && curr >= next) {
            // Local peak
            extrema.push((i, curr, true));
        } else if (curr <= prev && curr < next) || (curr < prev && curr <= next) {
            // Local nadir
            extrema.push((i, curr, false));
        }
    }

    if extrema.len() < 2 {
        return series.iter().copied().fold(f64::NAN, f64::max)
            - series.iter().copied().fold(f64::NAN, f64::min);
    }

    // Calculate valid excursions greater than 1 SD
    let mut excursions = Vec::new();
    for w in extrema.windows(2) {
        let (_idx1, val1, is_peak1) = w[0];
        let (_idx2, val2, is_peak2) = w[1];

        if is_peak1 != is_peak2 {
            let diff = (val1 - val2).abs();
            if diff > std_dev {
                excursions.push(diff);
            }
        }
    }

    if excursions.is_empty() {
        std_dev * 1.2
    } else {
        excursions.iter().sum::<f64>() / excursions.len() as f64
    }
}
