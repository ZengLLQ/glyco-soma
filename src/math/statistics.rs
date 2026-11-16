//! Statistical metrics, percentiles, distribution estimators, and sensor error evaluations.

use crate::error::{GlycoError, Result};

/// Comprehensive summary statistics for a biological time series.
#[derive(Debug, Clone, Default)]
pub struct SeriesStats {
    /// Number of observations.
    pub count: usize,
    /// Arithmetic mean.
    pub mean: f64,
    /// Sample standard deviation.
    pub std_dev: f64,
    /// Minimum value.
    pub min: f64,
    /// 5th percentile.
    pub p5: f64,
    /// 25th percentile (first quartile Q1).
    pub p25: f64,
    /// Median (50th percentile Q2).
    pub median: f64,
    /// 75th percentile (third quartile Q3).
    pub p75: f64,
    /// 95th percentile.
    pub p95: f64,
    /// Maximum value.
    pub max: f64,
    /// Interquartile Range (IQR = Q3 - Q1).
    pub iqr: f64,
    /// Coefficient of variation (CV% = std_dev / mean * 100).
    pub cv_percent: f64,
}

/// Calculate comprehensive summary statistics from a slice of floating point numbers.
pub fn compute_series_stats(data: &[f64]) -> Result<SeriesStats> {
    if data.is_empty() {
        return Err(GlycoError::EmptyDataSeries {
            series_name: "unnamed".to_string(),
        });
    }

    let n = data.len();
    let sum: f64 = data.iter().sum();
    let mean = sum / n as f64;

    let variance = if n > 1 {
        data.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1) as f64
    } else {
        0.0
    };
    let std_dev = variance.sqrt();

    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let min = sorted[0];
    let max = sorted[n - 1];

    let p5 = percentile_sorted(&sorted, 0.05);
    let p25 = percentile_sorted(&sorted, 0.25);
    let median = percentile_sorted(&sorted, 0.50);
    let p75 = percentile_sorted(&sorted, 0.75);
    let p95 = percentile_sorted(&sorted, 0.95);
    let iqr = p75 - p25;
    let cv_percent = if mean.abs() > 1e-9 {
        (std_dev / mean) * 100.0
    } else {
        0.0
    };

    Ok(SeriesStats {
        count: n,
        mean,
        std_dev,
        min,
        p5,
        p25,
        median,
        p75,
        p95,
        max,
        iqr,
        cv_percent,
    })
}

/// Linear interpolation percentile helper on pre-sorted array.
pub fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return sorted[0];
    }

    let rank = p.clamp(0.0, 1.0) * (n - 1) as f64;
    let low = rank.floor() as usize;
    let high = rank.ceil() as usize;
    let frac = rank - low as f64;

    if low == high {
        sorted[low]
    } else {
        sorted[low] * (1.0 - frac) + sorted[high] * frac
    }
}

/// Compute Mean Absolute Relative Difference (MARD) between sensor measurements and reference values:
///
/// $$ MARD = \frac{1}{N} \sum_{i=1}^N \frac{|y_{meas, i} - y_{ref, i}|}{y_{ref, i}} \times 100\% $$
pub fn compute_mard(meas: &[f64], reference: &[f64]) -> Result<f64> {
    if meas.len() != reference.len() {
        return Err(GlycoError::InvalidParameter {
            parameter: "reference.len".to_string(),
            value: reference.len() as f64,
            reason: format!(
                "Mismatched lengths: meas {}, ref {}",
                meas.len(),
                reference.len()
            ),
        });
    }
    if meas.is_empty() {
        return Err(GlycoError::EmptyDataSeries {
            series_name: "cgm_comparison".to_string(),
        });
    }

    let mut sum_ard = 0.0;
    let mut valid_count = 0;

    for (&m, &r) in meas.iter().zip(reference.iter()) {
        if r > 1e-6 {
            sum_ard += ((m - r).abs() / r) * 100.0;
            valid_count += 1;
        }
    }

    if valid_count == 0 {
        return Ok(0.0);
    }
    Ok(sum_ard / valid_count as f64)
}

/// Root-Mean-Square Error (RMSE) between two signals.
pub fn compute_rmse(a: &[f64], b: &[f64]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(GlycoError::InvalidParameter {
            parameter: "lengths".to_string(),
            value: a.len() as f64,
            reason: format!("Lengths do not match: {} vs {}", a.len(), b.len()),
        });
    }
    if a.is_empty() {
        return Err(GlycoError::EmptyDataSeries {
            series_name: "rmse".to_string(),
        });
    }

    let sum_sq: f64 = a
        .iter()
        .zip(b.iter())
        .map(|(&x, &y)| (x - y) * (x - y))
        .sum();
    Ok((sum_sq / a.len() as f64).sqrt())
}

/// Pearson product-moment correlation coefficient $r \in [-1.0, 1.0]$.
pub fn compute_pearson_r(x: &[f64], y: &[f64]) -> Result<f64> {
    if x.len() != y.len() || x.is_empty() {
        return Err(GlycoError::InvalidParameter {
            parameter: "vectors".to_string(),
            value: x.len() as f64,
            reason: "Vector lengths invalid or empty".to_string(),
        });
    }

    let n = x.len() as f64;
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;

    let mut cov = 0.0;
    let mut var_x = 0.0;
    let mut var_y = 0.0;

    for (&xi, &yi) in x.iter().zip(y.iter()) {
        let dx = xi - mean_x;
        let dy = yi - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    let denom = (var_x * var_y).sqrt();
    if denom < 1e-12 {
        return Ok(0.0);
    }
    Ok((cov / denom).clamp(-1.0, 1.0))
}

/// Gaussian Kernel Density Estimator for probability density evaluation.
#[derive(Debug, Clone)]
pub struct GaussianKde {
    samples: Vec<f64>,
    bandwidth: f64,
}

impl GaussianKde {
    /// Fit Gaussian KDE using Silverman's rule of thumb for optimal bandwidth.
    pub fn fit(samples: &[f64]) -> Result<Self> {
        if samples.len() < 2 {
            return Err(GlycoError::InvalidParameter {
                parameter: "samples.len".to_string(),
                value: samples.len() as f64,
                reason: "At least 2 samples needed for KDE estimation".to_string(),
            });
        }

        let stats = compute_series_stats(samples)?;
        // Silverman's rule of thumb: h = 0.9 * min(std, iqr/1.34) * n^(-1/5)
        let n = samples.len() as f64;
        let spread = stats.std_dev.min(stats.iqr / 1.34).max(1e-4);
        let bandwidth = (0.9 * spread * n.powf(-0.2)).max(1e-3);

        Ok(Self {
            samples: samples.to_vec(),
            bandwidth,
        })
    }

    /// Evaluate estimated probability density at point `x`.
    pub fn evaluate(&self, x: f64) -> f64 {
        let n = self.samples.len() as f64;
        let h = self.bandwidth;
        let inv_sqrt_2pi = 0.3989422804014327;

        let mut sum_kernel = 0.0;
        for &xi in &self.samples {
            let u = (x - xi) / h;
            sum_kernel += inv_sqrt_2pi * (-0.5 * u * u).exp();
        }

        sum_kernel / (n * h)
    }
}
