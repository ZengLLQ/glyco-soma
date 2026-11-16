//! Pure 7-bit ASCII Clark Error Grid scatter plot and zone distribution visualizer.

use crate::sensors::clark_grid::{classify_clark_point, ClarkGridReport, ClarkZone};

/// Render a 2D ASCII Clarke Error Grid scatter plot.
///
/// # Arguments
/// * `ref_vals`: Array of reference blood glucose values (mg/dL).
/// * `meas_vals`: Array of CGM sensor measured glucose values (mg/dL).
/// * `report`: Pre-computed Clarke Error Grid statistical report.
pub fn render_ascii_clark_plot(
    ref_vals: &[f64],
    meas_vals: &[f64],
    report: &ClarkGridReport,
) -> String {
    let mut out = String::new();

    out.push_str("+=========================================================================+\n");
    out.push_str("|            CLARKE ERROR GRID ANALYSIS (CEGA) SCATTER PLOT               |\n");
    out.push_str("+=========================================================================+\n");
    out.push_str(&format!(
        "| Total Pairs: {:5} | Zone A: {:5.1}% | Zone B: {:5.1}% | A+B: {:5.1}% (>=95%) |\n",
        report.total_pairs, report.percent_a, report.percent_b, report.percent_a_plus_b
    ));
    out.push_str(&format!(
        "| Zone C: {:5.1}%   | Zone D: {:5.1}% | Zone E: {:5.1}% | ISO Status: {}      |\n",
        report.percent_c,
        report.percent_d,
        report.percent_e,
        if report.iso_compliant { "PASS" } else { "FAIL" }
    ));
    out.push_str("+-------------------------------------------------------------------------+\n");

    // Grid dimension: 40 columns wide, 20 rows high (glucose range 0 to 400 mg/dL)
    let grid_w = 40;
    let grid_h = 20;
    let max_g = 400.0;

    let mut grid = vec![vec![' '; grid_w]; grid_h];

    // Draw identity diagonal line y = x using '.'
    for c in 0..grid_w {
        let r = (grid_h - 1) - (c * (grid_h - 1)) / (grid_w - 1);
        grid[r][c] = '.';
    }

    // Scatter points on the grid
    let n = ref_vals.len().min(meas_vals.len());
    for i in 0..n {
        let r = ref_vals[i].clamp(0.0, max_g);
        let m = meas_vals[i].clamp(0.0, max_g);

        let col = ((r / max_g) * (grid_w - 1) as f64).round() as usize;
        let row_idx = (grid_h - 1) - ((m / max_g) * (grid_h - 1) as f64).round() as usize;

        let zone = classify_clark_point(r, m);
        let symbol = match zone {
            ClarkZone::ZoneA => 'A',
            ClarkZone::ZoneB => 'B',
            ClarkZone::ZoneC => 'C',
            ClarkZone::ZoneD => 'D',
            ClarkZone::ZoneE => 'E',
        };

        grid[row_idx.min(grid_h - 1)][col.min(grid_w - 1)] = symbol;
    }

    out.push_str("  Measured (mg/dL)\n");
    for r in 0..grid_h {
        let y_val = max_g - (r as f64 / (grid_h - 1) as f64) * max_g;
        let line: String = grid[r].iter().collect();
        out.push_str(&format!("{:3.0} |{}|\n", y_val, line));
    }

    out.push_str("      +----------------------------------------+\n");
    out.push_str("      0       100       200       300       400 Reference (mg/dL)\n");
    out.push_str("+-------------------------------------------------------------------------+\n");
    out.push_str(
        "  Zones: [A] Clinically Accurate  [B] Benign  [C] Overcorrect  [D/E] Dangerous\n",
    );
    out.push_str("+=========================================================================+\n");

    out
}
