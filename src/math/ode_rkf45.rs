//! Adaptive Runge-Kutta-Fehlberg (RKF45) numerical ODE integrator with local error control.

use crate::error::{GlycoError, Result};

/// Fehlberg 4(5) Butcher tableau coefficients.
const C: [f64; 6] = [0.0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0];

const A: [[f64; 5]; 5] = [
    [1.0 / 4.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0],
    [1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0],
    [439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0.0],
    [
        -8.0 / 27.0,
        2.0,
        -3544.0 / 2565.0,
        1859.0 / 4104.0,
        -11.0 / 40.0,
    ],
];

/// 4th-order weights.
const B4: [f64; 6] = [
    25.0 / 216.0,
    0.0,
    1408.0 / 2565.0,
    2197.0 / 4104.0,
    -1.0 / 5.0,
    0.0,
];

/// 5th-order weights.
const B5: [f64; 6] = [
    16.0 / 135.0,
    0.0,
    6656.0 / 12825.0,
    28561.0 / 56430.0,
    -9.0 / 50.0,
    2.0 / 55.0,
];

/// Configuration options for the RKF45 adaptive solver.
#[derive(Debug, Clone)]
pub struct Rkf45Config {
    /// Relative and absolute error tolerance $\epsilon$.
    pub tolerance: f64,
    /// Minimum permissible time step (minutes).
    pub min_dt: f64,
    /// Maximum permissible time step (minutes).
    pub max_dt: f64,
    /// Safety factor for step size adaptation (typically 0.85 to 0.90).
    pub safety_factor: f64,
    /// Maximum step size scaling ratio per step.
    pub max_scale: f64,
    /// Minimum step size scaling ratio per step.
    pub min_scale: f64,
}

impl Default for Rkf45Config {
    fn default() -> Self {
        Self {
            tolerance: 1e-5,
            min_dt: 1e-4,
            max_dt: 2.0,
            safety_factor: 0.88,
            max_scale: 2.0,
            min_scale: 0.2,
        }
    }
}

/// Statistics collected across an integration trajectory.
#[derive(Debug, Clone, Default)]
pub struct SolverStats {
    /// Total successful time steps taken.
    pub accepted_steps: usize,
    /// Steps rejected due to excessive local truncation error.
    pub rejected_steps: usize,
    /// Total function evaluations (f(t, y)).
    pub function_evals: usize,
}

/// Adaptive Runge-Kutta-Fehlberg 4(5) solver.
pub struct Rkf45Solver {
    config: Rkf45Config,
    pub stats: SolverStats,
}

impl Rkf45Solver {
    /// Construct a new RKF45 solver instance.
    pub fn new(config: Rkf45Config) -> Self {
        Self {
            config,
            stats: SolverStats::default(),
        }
    }

    /// Advance the state vector `y` from time `t` across interval `dt_target` using adaptive sub-stepping.
    ///
    /// The derivative function `f(t, y, dydt)` evaluates the vector field in place.
    pub fn integrate_step<F>(
        &mut self,
        t_start: f64,
        t_end: f64,
        y: &mut [f64],
        mut f: F,
    ) -> Result<f64>
    where
        F: FnMut(f64, &[f64], &mut [f64]),
    {
        let dim = y.len();
        let mut t = t_start;
        let mut h = (t_end - t_start).min(self.config.max_dt);

        let mut k = vec![vec![0.0; dim]; 6];
        let mut y_temp = vec![0.0; dim];
        let mut y4 = vec![0.0; dim];
        let mut y5 = vec![0.0; dim];

        while t < t_end - 1e-12 {
            // Clamp h so we land exactly at t_end
            if t + h > t_end {
                h = t_end - t;
            }

            // Stage 1
            f(t, y, &mut k[0]);
            self.stats.function_evals += 1;

            // Stages 2 to 6
            for s in 1..6 {
                for i in 0..dim {
                    let mut sum_ak = 0.0;
                    for j in 0..s {
                        sum_ak += A[s - 1][j] * k[j][i];
                    }
                    y_temp[i] = y[i] + h * sum_ak;
                }
                f(t + C[s] * h, &y_temp, &mut k[s]);
                self.stats.function_evals += 1;
            }

            // Compute 4th and 5th order approximations
            let mut err_max = 0.0;
            for i in 0..dim {
                let mut sum_b4 = 0.0;
                let mut sum_b5 = 0.0;
                for s in 0..6 {
                    sum_b4 += B4[s] * k[s][i];
                    sum_b5 += B5[s] * k[s][i];
                }
                y4[i] = y[i] + h * sum_b4;
                y5[i] = y[i] + h * sum_b5;

                // Scale error relative to magnitude
                let sc = 1.0 + y[i].abs().max(y5[i].abs());
                let err_i = (y5[i] - y4[i]).abs() / sc;
                if err_i > err_max {
                    err_max = err_i;
                }
            }

            let normalized_err = err_max / self.config.tolerance;

            if normalized_err <= 1.0 || h <= self.config.min_dt {
                // Step accepted
                t += h;
                for i in 0..dim {
                    // Positivity invariant enforcement for biological concentrations
                    y[i] = y5[i].max(0.0);
                }
                self.stats.accepted_steps += 1;

                // Calculate next optimal step size
                let scale = if normalized_err < 1e-12 {
                    self.config.max_scale
                } else {
                    (self.config.safety_factor * (1.0 / normalized_err).powf(0.2))
                        .clamp(self.config.min_scale, self.config.max_scale)
                };
                h = (h * scale).clamp(self.config.min_dt, self.config.max_dt);
            } else {
                // Step rejected, shrink step size
                self.stats.rejected_steps += 1;
                let scale = (self.config.safety_factor * (1.0 / normalized_err).powf(0.25))
                    .clamp(self.config.min_scale, 0.9);
                h *= scale;

                if h < self.config.min_dt {
                    if (t_end - t).abs() < self.config.min_dt * 2.0 {
                        // Very close to endpoint, force proceed
                        t = t_end;
                        for i in 0..dim {
                            y[i] = y4[i].max(0.0);
                        }
                        break;
                    } else {
                        return Err(GlycoError::TimeStepTooSmall {
                            dt: h,
                            min_dt: self.config.min_dt,
                            time: t,
                        });
                    }
                }
            }
        }

        Ok(t)
    }
}
