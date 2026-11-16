//! Pure 7-bit ASCII line chart and time series strip plotter for terminal displays.
//!
//! Strictly uses ASCII codes 0-127 (`*`, `-`, `+`, `|`, `#`, `.`, `v`, `^`, space)
//! so Windows cp1252 and headless CI environments never crash with `UnicodeEncodeError`.

/// Render a 2D time-series glucose chart using pure 7-bit ASCII characters.
///
/// # Arguments
/// * `time_min`: Array of time points (minutes).
/// * `glucose_mg_dl`: Array of glucose values (mg/dL).
/// * `width`: Terminal column width (e.g. 70 to 90).
/// * `height`: Terminal row height (e.g. 18 to 25).
/// * `target_low`: Lower target bound (typically 70.0 mg/dL).
/// * `target_high`: Upper target bound (typically 180.0 mg/dL).
pub fn render_ascii_glucose_chart(
    time_min: &[f64],
    glucose_mg_dl: &[f64],
    width: usize,
    height: usize,
    target_low: f64,
    target_high: f64,
) -> String {
    if glucose_mg_dl.is_empty() || time_min.is_empty() {
        return "No data available to plot.".to_string();
    }

    let n = glucose_mg_dl.len().min(time_min.len());
    let plot_w = width.clamp(30, 120);
    let plot_h = height.clamp(10, 40);

    // Dynamic vertical bounds
    let data_min = glucose_mg_dl.iter().copied().fold(f64::INFINITY, f64::min);
    let data_max = glucose_mg_dl
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    let y_min = (data_min.min(target_low) - 10.0).max(30.0).floor();
    let y_max = (data_max.max(target_high) + 15.0).min(450.0).ceil();

    let t_min = time_min[0];
    let t_max = time_min[n - 1].max(t_min + 1.0);

    // Initialize 2D grid with spaces
    let mut grid = vec![vec![' '; plot_w]; plot_h];

    // Compute row indices for target thresholds (row 0 is top / y_max, row plot_h - 1 is bottom / y_min)
    let y_to_row = |y: f64| -> usize {
        let frac = (y_max - y) / (y_max - y_min).max(1.0);
        let row = (frac * (plot_h - 1) as f64).round() as isize;
        row.clamp(0, (plot_h - 1) as isize) as usize
    };

    let row_target_high = y_to_row(target_high);
    let row_target_low = y_to_row(target_low);

    // Draw horizontal target threshold lines using '-'
    for x in 0..plot_w {
        if grid[row_target_high][x] == ' ' {
            grid[row_target_high][x] = '-';
        }
        if grid[row_target_low][x] == ' ' {
            grid[row_target_low][x] = '-';
        }
    }

    // Plot data points using '*'
    for i in 0..n {
        let t = time_min[i];
        let g = glucose_mg_dl[i];

        let col_frac = (t - t_min) / (t_max - t_min).max(1.0);
        let col = ((col_frac * (plot_w - 1) as f64).round() as usize).min(plot_w - 1);
        let row = y_to_row(g);

        // Character based on clinical severity
        let ch = if g < 54.0 {
            '!' // Urgent hypoglycemia
        } else if g < 70.0 {
            'v' // Hypoglycemia
        } else if g > 250.0 {
            '#' // Severe hyperglycemia
        } else if g > 180.0 {
            '^' // Hyperglycemia
        } else {
            '*' // In target range
        };

        grid[row][col] = ch;
    }

    // Format output with y-axis labels and left border '|'
    let mut out = String::new();
    out.push_str("+-------------------------------------------------------------------------+\n");
    out.push_str("| GLUCOSE STRIP CHART (70-180 mg/dL target indicated by '-')              |\n");
    out.push_str("+-------------------------------------------------------------------------+\n");

    for r in 0..plot_h {
        let y_val = y_max - (r as f64 / (plot_h - 1) as f64) * (y_max - y_min);
        let label = if r == 0 || r == plot_h - 1 || r == row_target_high || r == row_target_low {
            format!("{:3.0}", y_val)
        } else {
            "   ".to_string()
        };

        let row_str: String = grid[r].iter().collect();
        out.push_str(&format!("{} |{}|\n", label, row_str));
    }

    // Bottom horizontal border and x-axis labels
    let border: String = std::iter::repeat('-').take(plot_w).collect();
    out.push_str(&format!("    +{}\n", border));

    // Time markings
    let label_t0 = format!("{:.0}m", t_min);
    let label_tmid = format!("{:.0}m", (t_min + t_max) / 2.0);
    let label_tend = format!("{:.0}m", t_max);

    let space_mid = (plot_w / 2).saturating_sub(label_t0.len() + label_tmid.len() / 2);
    let space_end =
        plot_w.saturating_sub(label_t0.len() + space_mid + label_tmid.len() + label_tend.len());

    let pad_mid: String = std::iter::repeat(' ').take(space_mid).collect();
    let pad_end: String = std::iter::repeat(' ').take(space_end).collect();

    out.push_str(&format!(
        "     {}{}{}{}\n",
        label_t0, pad_mid, label_tmid, pad_end
    ));
    out.push_str("     Legend: [*] Target  [v] Hypo (<70)  [!] Severe (<54)  [^] Hyper (>180)  [#] High (>250)\n");

    out
}
