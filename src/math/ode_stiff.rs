//! Implicit second-order Crank-Nicolson and Backward Euler solvers with Newton-Raphson iterations for stiff metabolic dynamics.

use crate::error::{GlycoError, Result};

/// Stiff ODE solver configuration.
#[derive(Debug, Clone)]
pub struct StiffSolverConfig {
    /// Maximum Newton-Raphson non-linear iterations per time step.
    pub max_newton_iters: usize,
    /// Absolute convergence tolerance on the residual norm.
    pub tolerance: f64,
    /// Finite-difference perturbation magnitude for numerical Jacobian calculation.
    pub fd_eps: f64,
    /// Crank-Nicolson theta weighting (0.5 for second-order A-stable Crank-Nicolson, 1.0 for L-stable Backward Euler).
    pub theta: f64,
}

impl Default for StiffSolverConfig {
    fn default() -> Self {
        Self {
            max_newton_iters: 25,
            tolerance: 1e-6,
            fd_eps: 1e-7,
            theta: 0.5, // Crank-Nicolson
        }
    }
}

/// Stiff implicit differential equation solver.
pub struct StiffSolver {
    config: StiffSolverConfig,
}

impl StiffSolver {
    /// Create a new solver instance with given configuration.
    pub fn new(config: StiffSolverConfig) -> Self {
        Self { config }
    }

    /// Advance state `y` from time `t` to `t + dt` using theta-method implicit integration:
    ///
    /// $$ y_{n+1} - y_n - dt \left[ (1-\theta) f(t_n, y_n) + \theta f(t_{n+1}, y_{n+1}) \right] = 0 $$
    pub fn step<F>(&self, t: f64, dt: f64, y: &mut [f64], mut f: F) -> Result<()>
    where
        F: FnMut(f64, &[f64], &mut [f64]),
    {
        let dim = y.len();
        let mut fn_val = vec![0.0; dim];
        f(t, y, &mut fn_val);

        let t_next = t + dt;
        let theta = self.config.theta;
        let mut y_next = y.to_vec(); // Initial guess: y_{n+1}^{(0)} = y_n

        let mut res = vec![0.0; dim];
        let mut f_next = vec![0.0; dim];
        let mut jac = vec![vec![0.0; dim]; dim];
        let mut delta_y = vec![0.0; dim];

        for iter in 0..self.config.max_newton_iters {
            f(t_next, &y_next, &mut f_next);

            // Compute residual: R = y_next - y - dt * ((1-theta)*fn + theta*f_next)
            let mut res_norm = 0.0;
            for i in 0..dim {
                res[i] = y_next[i] - y[i] - dt * ((1.0 - theta) * fn_val[i] + theta * f_next[i]);
                res_norm += res[i] * res[i];
            }
            res_norm = res_norm.sqrt();

            if res_norm < self.config.tolerance {
                // Converged successfully
                for i in 0..dim {
                    y[i] = y_next[i].max(0.0);
                }
                return Ok(());
            }

            // Assemble Jacobian: J = I - dt * theta * \partial f / \partial y
            for j in 0..dim {
                let y_orig = y_next[j];
                let h_fd = self.config.fd_eps * (1.0 + y_orig.abs());
                y_next[j] += h_fd;

                let mut f_pert = vec![0.0; dim];
                f(t_next, &y_next, &mut f_pert);
                y_next[j] = y_orig;

                for i in 0..dim {
                    let df_dy = (f_pert[i] - f_next[i]) / h_fd;
                    let identity = if i == j { 1.0 } else { 0.0 };
                    jac[i][j] = identity - dt * theta * df_dy;
                }
            }

            // Solve linear system: J * delta_y = -res via Gaussian elimination with partial pivoting
            solve_linear_system(&mut jac, &res, &mut delta_y)?;

            // Update state with damping line search
            let mut alpha = 1.0;
            let mut improved = false;
            for _ in 0..5 {
                let mut candidate_y = vec![0.0; dim];
                for i in 0..dim {
                    candidate_y[i] = (y_next[i] - alpha * delta_y[i]).max(0.0);
                }

                let mut f_cand = vec![0.0; dim];
                f(t_next, &candidate_y, &mut f_cand);

                let mut cand_res_norm = 0.0;
                for i in 0..dim {
                    let r = candidate_y[i]
                        - y[i]
                        - dt * ((1.0 - theta) * fn_val[i] + theta * f_cand[i]);
                    cand_res_norm += r * r;
                }
                cand_res_norm = cand_res_norm.sqrt();

                if cand_res_norm < res_norm {
                    y_next = candidate_y;
                    improved = true;
                    break;
                }
                alpha *= 0.5;
            }

            if !improved {
                for i in 0..dim {
                    y_next[i] = (y_next[i] - 0.25 * delta_y[i]).max(0.0);
                }
            }

            if iter == self.config.max_newton_iters - 1 {
                return Err(GlycoError::OptimizationFailed {
                    iterations: iter,
                    residual: res_norm,
                    reason: "Newton-Raphson iteration exceeded maximum steps in stiff solver"
                        .to_string(),
                });
            }
        }

        for i in 0..dim {
            y[i] = y_next[i].max(0.0);
        }
        Ok(())
    }
}

/// Gaussian elimination solver with partial pivoting for moderate dimensions ($N \le 20$).
fn solve_linear_system(a: &mut [Vec<f64>], b: &[f64], x: &mut [f64]) -> Result<()> {
    let n = b.len();
    let mut rhs = b.to_vec();

    for i in 0..n {
        // Find pivot
        let mut max_row = i;
        let mut max_val = a[i][i].abs();
        for k in (i + 1)..n {
            if a[k][i].abs() > max_val {
                max_val = a[k][i].abs();
                max_row = k;
            }
        }

        if max_val < 1e-13 {
            return Err(GlycoError::MatrixSingular {
                dimension: n,
                detail: format!("Zero pivot encountered at row {}", i),
            });
        }

        // Swap rows
        if max_row != i {
            a.swap(i, max_row);
            rhs.swap(i, max_row);
        }

        // Eliminate
        for k in (i + 1)..n {
            let factor = a[k][i] / a[i][i];
            for j in (i + 1)..n {
                a[k][j] -= factor * a[i][j];
            }
            rhs[k] -= factor * rhs[i];
            a[k][i] = 0.0;
        }
    }

    // Back-substitution
    for i in (0..n).rev() {
        let mut sum = rhs[i];
        for j in (i + 1)..n {
            sum -= a[i][j] * x[j];
        }
        x[i] = sum / a[i][i];
    }

    Ok(())
}
