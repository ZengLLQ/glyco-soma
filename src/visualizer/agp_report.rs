//! Ambulatory Glucose Profile (AGP) 24-hour standardized clinical report in pure 7-bit ASCII.

use crate::clinical::attd_metrics::AttdMetricsReport;
use crate::clinical::variability::GlycemicVariabilityReport;

/// Render a 24-hour Ambulatory Glucose Profile (AGP) report in pure 7-bit ASCII.
pub fn render_ascii_agp_report(
    attd: &AttdMetricsReport,
    variability: &GlycemicVariabilityReport,
    hourly_medians: &[f64; 24],
    hourly_p25: &[f64; 24],
    hourly_p75: &[f64; 24],
) -> String {
    let mut s = String::new();

    s.push_str("+=========================================================================+\n");
    s.push_str("|            AMBULATORY GLUCOSE PROFILE (AGP) CLINICAL REPORT             |\n");
    s.push_str("+=========================================================================+\n");
    s.push_str(&format!(
        "| Time Analyzed: {:5.1} hours | Mean Glucose: {:5.1} mg/dL | GMI: {:4.1}% (HbA1c)    |\n",
        attd.total_hours, attd.mean_glucose_mg_dl, attd.gmi_percent
    ));
    s.push_str(&format!(
        "| Glucose CV:    {:5.1}%       | Std Dev:      {:5.1} mg/dL | MAGE: {:5.1} mg/dL    |\n",
        variability.cv_percent, variability.std_dev_mg_dl, variability.mage_mg_dl
    ));
    s.push_str("+-------------------------------------------------------------------------+\n");
    s.push_str("| TIME IN RANGES (ATTD 2024 International Consensus Target Goals)         |\n");
    s.push_str("+-------------------------------------------------------------------------+\n");
    s.push_str(&format!(
        "| Very High (>250 mg/dL, TAR L2): {:5.1}%  (Target < 5%)   {}\n",
        attd.tar_l2_percent,
        render_bar(attd.tar_l2_percent, 100.0, 20, '#')
    ));
    s.push_str(&format!(
        "| High (181-250 mg/dL,   TAR L1): {:5.1}%  (Target < 25%)  {}\n",
        attd.tar_l1_percent,
        render_bar(attd.tar_l1_percent, 100.0, 20, '^')
    ));
    s.push_str(&format!(
        "| In Target (70-180 mg/dL,  TIR): {:5.1}%  (Target > 70%)  {}\n",
        attd.tir_percent,
        render_bar(attd.tir_percent, 100.0, 20, '*')
    ));
    s.push_str(&format!(
        "| Low (54-69 mg/dL,      TBR L1): {:5.1}%  (Target < 4%)   {}\n",
        attd.tbr_l1_percent,
        render_bar(attd.tbr_l1_percent, 100.0, 20, 'v')
    ));
    s.push_str(&format!(
        "| Very Low (<54 mg/dL,   TBR L2): {:5.1}%  (Target < 1%)   {}\n",
        attd.tbr_l2_percent,
        render_bar(attd.tbr_l2_percent, 100.0, 20, '!')
    ));
    s.push_str("+-------------------------------------------------------------------------+\n");
    s.push_str(&format!(
        "| Glycemic Risk Index (GRI): {:4.1}/100 [{}] |\n",
        attd.glycemic_risk_index, attd.gri_category
    ));
    s.push_str(&format!(
        "| Hypoglycemia Index (LBGI): {:4.1} [{}]                        |\n",
        variability.lbgi, variability.lbgi_risk_category
    ));
    s.push_str(&format!(
        "| Hyperglycemia Index (HBGI):{:4.1} [{}]                        |\n",
        variability.hbgi, variability.hbgi_risk_category
    ));
    s.push_str("+-------------------------------------------------------------------------+\n");
    s.push_str("| 24-HOUR MODAL DAY PROFILE (Median [*] & 25th-75th Percentile Band [:])  |\n");
    s.push_str("+-------------------------------------------------------------------------+\n");

    // 24-hour grid: 24 columns (one per hour), 12 rows (vertical glucose scale 50 to 250)
    let grid_h = 12;
    let grid_w = 24;
    let y_min = 50.0;
    let y_max = 250.0;

    let mut grid = vec![vec![' '; grid_w * 2]; grid_h];

    let row_at = |val: f64| -> usize {
        let frac = (y_max - val) / (y_max - y_min);
        let r = (frac * (grid_h - 1) as f64).round() as isize;
        r.clamp(0, (grid_h - 1) as isize) as usize
    };

    let row_180 = row_at(180.0);
    let row_70 = row_at(70.0);

    // Target threshold lines
    for x in 0..(grid_w * 2) {
        grid[row_180][x] = '-';
        grid[row_70][x] = '-';
    }

    for h in 0..24 {
        let col = h * 2;
        let r_med = row_at(hourly_medians[h]);
        let r_p25 = row_at(hourly_p25[h]);
        let r_p75 = row_at(hourly_p75[h]);

        // Draw interquartile band
        let top_r = r_p75.min(r_p25);
        let bot_r = r_p75.max(r_p25);
        for r in top_r..=bot_r {
            if grid[r][col] == ' ' || grid[r][col] == '-' {
                grid[r][col] = ':';
            }
        }

        // Draw median point
        grid[r_med][col] = '*';
    }

    for r in 0..grid_h {
        let y_val = y_max - (r as f64 / (grid_h - 1) as f64) * (y_max - y_min);
        let line: String = grid[r].iter().collect();
        s.push_str(&format!("{:3.0} |{}|\n", y_val, line));
    }

    s.push_str("    +------------------------------------------------+\n");
    s.push_str("     00  02  04  06  08  10  12  14  16  18  20  22 Hour\n");
    s.push_str("+=========================================================================+\n");

    s
}

fn render_bar(val: f64, max: f64, bar_len: usize, ch: char) -> String {
    let filled = ((val / max) * bar_len as f64).round() as usize;
    let filled_clamped = filled.min(bar_len);
    let empty = bar_len - filled_clamped;
    let f_str: String = std::iter::repeat(ch).take(filled_clamped).collect();
    let e_str: String = std::iter::repeat('.').take(empty).collect();
    format!("[{}{}]", f_str, e_str)
}
