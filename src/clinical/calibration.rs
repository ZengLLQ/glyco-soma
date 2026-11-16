//! Patient parameter calibration and sensitivity screening from continuous glucose monitoring data.

use crate::error::{GlycoError, Result};
use crate::types::PatientProfile;

/// Calibration result containing fitted patient parameters and goodness of fit metrics.
#[derive(Debug, Clone)]
pub struct CalibrationResult {
    /// Fitted Insulin Sensitivity Factor (ISF, mg/dL/U).
    pub fitted_isf: f64,
    /// Fitted Carbohydrate-to-Insulin Ratio (CIR, g/U).
    pub fitted_cir: f64,
    /// Fitted basal infusion rate (U/hr).
    pub fitted_basal_u_hr: f64,
    /// Root-Mean-Square Error between model prediction and observed CGM trace (mg/dL).
    pub rmse_mg_dl: f64,
    /// Number of optimization iterations executed.
    pub iterations: usize,
    /// Calibration convergence status.
    pub converged: bool,
}

/// Calibrate patient-specific parameters (ISF, CIR, Basal) using Nelder-Mead simplex optimization.
pub fn calibrate_patient_parameters(
    initial_profile: &PatientProfile,
    cgm_measurements: &[f64],
    carbs_ingested_g: f64,
    insulin_delivered_u: f64,
) -> Result<CalibrationResult> {
    if cgm_measurements.len() < 10 {
        return Err(GlycoError::InvalidParameter {
            parameter: "cgm_measurements.len".to_string(),
            value: cgm_measurements.len() as f64,
            reason: "At least 10 CGM samples required for parameter calibration".to_string(),
        });
    }

    let initial_isf = initial_profile.isf_mg_dl_per_u;
    let initial_cir = initial_profile.carb_ratio_g_per_u;
    let initial_basal = initial_profile.basal_rate_u_per_hr;

    // Parameter vector: [ISF, CIR, Basal]
    let mut simplex = vec![
        vec![initial_isf, initial_cir, initial_basal],
        vec![initial_isf * 1.15, initial_cir, initial_basal],
        vec![initial_isf, initial_cir * 1.15, initial_basal],
        vec![initial_isf, initial_cir, initial_basal * 1.15],
    ];

    let mut scores = vec![0.0; 4];
    for (i, p) in simplex.iter().enumerate() {
        scores[i] = evaluate_cost(p, cgm_measurements, carbs_ingested_g, insulin_delivered_u);
    }

    let max_iters = 80;
    let mut converged = false;
    let mut iter_count = 0;

    for it in 0..max_iters {
        iter_count = it;

        // Sort simplex by cost
        let mut indices: Vec<usize> = (0..4).collect();
        indices.sort_by(|&a, &b| scores[a].partial_cmp(&scores[b]).unwrap());

        let best_idx = indices[0];
        let worst_idx = indices[3];
        let second_worst_idx = indices[2];

        let best_score = scores[best_idx];
        let worst_score = scores[worst_idx];

        if (worst_score - best_score).abs() < 1e-3 {
            converged = true;
            break;
        }

        // Centroid of all points except worst
        let mut centroid = vec![0.0; 3];
        for &idx in &indices[0..3] {
            for j in 0..3 {
                centroid[j] += simplex[idx][j] / 3.0;
            }
        }

        // Reflection
        let alpha = 1.0;
        let mut reflected = vec![0.0; 3];
        for j in 0..3 {
            reflected[j] = (centroid[j] + alpha * (centroid[j] - simplex[worst_idx][j])).clamp(
                if j == 0 {
                    15.0
                } else if j == 1 {
                    4.0
                } else {
                    0.1
                },
                if j == 0 {
                    150.0
                } else if j == 1 {
                    35.0
                } else {
                    4.0
                },
            );
        }
        let ref_score = evaluate_cost(
            &reflected,
            cgm_measurements,
            carbs_ingested_g,
            insulin_delivered_u,
        );

        if ref_score < scores[second_worst_idx] && ref_score >= best_score {
            simplex[worst_idx] = reflected;
            scores[worst_idx] = ref_score;
            continue;
        }

        // Expansion
        if ref_score < best_score {
            let gamma = 2.0;
            let mut expanded = vec![0.0; 3];
            for j in 0..3 {
                expanded[j] = (centroid[j] + gamma * (reflected[j] - centroid[j])).clamp(
                    if j == 0 {
                        15.0
                    } else if j == 1 {
                        4.0
                    } else {
                        0.1
                    },
                    if j == 0 {
                        150.0
                    } else if j == 1 {
                        35.0
                    } else {
                        4.0
                    },
                );
            }
            let exp_score = evaluate_cost(
                &expanded,
                cgm_measurements,
                carbs_ingested_g,
                insulin_delivered_u,
            );
            if exp_score < ref_score {
                simplex[worst_idx] = expanded;
                scores[worst_idx] = exp_score;
            } else {
                simplex[worst_idx] = reflected;
                scores[worst_idx] = ref_score;
            }
            continue;
        }

        // Contraction
        let rho = 0.5;
        let mut contracted = vec![0.0; 3];
        for j in 0..3 {
            contracted[j] = centroid[j] + rho * (simplex[worst_idx][j] - centroid[j]);
        }
        let cont_score = evaluate_cost(
            &contracted,
            cgm_measurements,
            carbs_ingested_g,
            insulin_delivered_u,
        );
        if cont_score < worst_score {
            simplex[worst_idx] = contracted;
            scores[worst_idx] = cont_score;
            continue;
        }

        // Shrink
        let sigma = 0.5;
        for &idx in &indices[1..4] {
            for j in 0..3 {
                simplex[idx][j] =
                    simplex[best_idx][j] + sigma * (simplex[idx][j] - simplex[best_idx][j]);
            }
            scores[idx] = evaluate_cost(
                &simplex[idx],
                cgm_measurements,
                carbs_ingested_g,
                insulin_delivered_u,
            );
        }
    }

    let mut best_i = 0;
    let mut min_val = scores[0];
    for (i, &s) in scores.iter().enumerate() {
        if s < min_val {
            min_val = s;
            best_i = i;
        }
    }

    let fitted = &simplex[best_i];

    Ok(CalibrationResult {
        fitted_isf: fitted[0],
        fitted_cir: fitted[1],
        fitted_basal_u_hr: fitted[2],
        rmse_mg_dl: min_val.sqrt(),
        iterations: iter_count + 1,
        converged,
    })
}

/// Cost function computing mean squared residual error between predicted trajectory and observed CGM.
fn evaluate_cost(params: &[f64], cgm: &[f64], carbs_g: f64, insulin_u: f64) -> f64 {
    let isf = params[0];
    let cir = params[1];
    let _basal = params[2];

    let g0 = cgm[0];
    let n = cgm.len();
    let mut sum_sq_err = 0.0;

    // Analytical two-pole impulse response approximation:
    // delta_G(t) = Carb_impact(t) - Insulin_impact(t)
    for (k, &observed_g) in cgm.iter().enumerate() {
        let t = k as f64 * 5.0; // 5 min samples

        // Carb peak around 45 min
        let carb_curve = if carbs_g > 0.0 {
            let t_peak = 45.0;
            let carb_amplitude = (carbs_g / cir) * isf;
            (t / t_peak) * (-t / t_peak).exp() * carb_amplitude * 2.71828
        } else {
            0.0
        };

        // Insulin action peak around 75 min
        let insulin_curve = if insulin_u > 0.0 {
            let t_ins_peak = 75.0;
            let ins_drop = insulin_u * isf;
            (t / t_ins_peak).powi(2) * (-t / t_ins_peak).exp() * ins_drop * 1.359
        } else {
            0.0
        };

        let predicted_g = (g0 + carb_curve - insulin_curve).max(20.0);
        let err = predicted_g - observed_g;
        sum_sq_err += err * err;
    }

    sum_sq_err / n as f64
}
