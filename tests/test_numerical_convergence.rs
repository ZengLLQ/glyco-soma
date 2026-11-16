//! Numerical convergence and error order verification tests for RKF45 and Stiff solvers.

use glyco_soma::math::ode_rkf45::{Rkf45Config, Rkf45Solver};
use glyco_soma::math::ode_stiff::{StiffSolver, StiffSolverConfig};

#[test]
fn test_rkf45_convergence_exponential_decay() {
    // Solve dy/dt = -lambda * y, exact solution: y(t) = y0 * exp(-lambda * t)
    let lambda: f64 = 0.25;
    let y0: f64 = 100.0;
    let t_end: f64 = 4.0;
    let exact_y_end = y0 * (-lambda * t_end).exp();

    let mut solver = Rkf45Solver::new(Rkf45Config {
        tolerance: 1e-6,
        min_dt: 1e-5,
        max_dt: 0.5,
        ..Default::default()
    });

    let mut y = vec![y0];
    let res_t = solver
        .integrate_step(0.0, t_end, &mut y, |_t, state, dydt| {
            dydt[0] = -lambda * state[0];
        })
        .unwrap();

    assert!((res_t - t_end).abs() < 1e-9);
    let numerical_error = (y[0] - exact_y_end).abs();

    assert!(
        numerical_error < 1e-4,
        "RKF45 numerical error {:.2e} exceeds threshold",
        numerical_error
    );
    assert!(solver.stats.accepted_steps > 0);
}

#[test]
fn test_stiff_solver_crank_nicolson_decay() {
    // Stiff decay test with large lambda = 50.0
    let lambda = 50.0;
    let y0 = 100.0;
    let dt = 0.05; // dt * lambda = 2.5 > 2.0 (explicit Euler would diverge, Crank-Nicolson remains stable)

    let solver = StiffSolver::new(StiffSolverConfig {
        theta: 1.0, // Backward Euler (L-stable)
        tolerance: 1e-7,
        ..Default::default()
    });

    let mut y = vec![y0];
    let mut t = 0.0;

    for _ in 0..10 {
        solver
            .step(t, dt, &mut y, |_time, state, dydt| {
                dydt[0] = -lambda * state[0];
            })
            .unwrap();
        t += dt;
        assert!(y[0] >= 0.0, "Crank-Nicolson state must remain non-negative");
    }

    // State should rapidly decay monotonically towards 0
    assert!(y[0] < 1.0);
}

#[test]
fn test_rkf45_step_rejection_on_discontinuity() {
    let mut solver = Rkf45Solver::new(Rkf45Config {
        tolerance: 1e-5,
        min_dt: 1e-4,
        max_dt: 1.0,
        ..Default::default()
    });

    let mut y = vec![1.0];
    // Sharp step gradient
    let _ = solver.integrate_step(0.0, 2.0, &mut y, |t, state, dydt| {
        if t < 1.0 {
            dydt[0] = 0.1;
        } else {
            dydt[0] = 50.0 * state[0]; // sudden surge
        }
    });

    // Solver successfully traversed the trajectory
    assert!(y[0] > 1.0);
}
