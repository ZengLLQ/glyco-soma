//! Pure 7-bit ASCII visualizers, clinical reports, and Clarke Error Grid scatter plots.

pub mod agp_report;
pub mod ascii_chart;
pub mod clark_plot;

pub use agp_report::render_ascii_agp_report;
pub use ascii_chart::render_ascii_glucose_chart;
pub use clark_plot::render_ascii_clark_plot;
