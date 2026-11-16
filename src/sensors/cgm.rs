//! Subcutaneous Continuous Glucose Monitor (CGM) simulator with diffusion delay, colored AR(1) noise, and compression artifacts.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Normal};

/// CGM hardware sensor and physiological parameters.
#[derive(Debug, Clone)]
pub struct CgmParams {
    /// Interstitial-to-vascular glucose diffusion time constant $\tau_{cgm}$ in minutes (typically 8.0 to 12.0 min).
    pub tau_interstitial_min: f64,
    /// Autoregressive AR(1) noise correlation coefficient $\phi \in [0.7, 0.95]$.
    pub ar1_phi: f64,
    /// White noise innovation standard deviation in mg/dL.
    pub noise_std_mg_dl: f64,
    /// Sensor linear calibration drift slope (% per day).
    pub drift_percent_per_day: f64,
    /// Probability of a compression artifact per hour (e.g. nocturnal sleeping on sensor).
    pub compression_low_prob_per_hr: f64,
    /// Maximum sensor reading saturation limit (mg/dL).
    pub max_reading_mg_dl: f64,
    /// Minimum sensor reading lower limit (mg/dL).
    pub min_reading_mg_dl: f64,
    /// Sampling period in minutes (typically 5.0 min).
    pub sampling_interval_min: f64,
}

impl Default for CgmParams {
    fn default() -> Self {
        Self {
            tau_interstitial_min: 9.5,
            ar1_phi: 0.82,
            noise_std_mg_dl: 3.5,
            drift_percent_per_day: 1.5,
            compression_low_prob_per_hr: 0.03,
            max_reading_mg_dl: 400.0,
            min_reading_mg_dl: 40.0,
            sampling_interval_min: 5.0,
        }
    }
}

/// Dynamic state of the CGM sensor.
#[derive(Debug, Clone)]
pub struct CgmState {
    /// Interstitial glucose concentration $G_{sub}(t)$ (mg/dL).
    pub interstitial_glucose_mg_dl: f64,
    /// Current AR(1) correlated noise component $\epsilon_k$.
    pub ar1_noise_val: f64,
    /// Active compression artifact remaining duration (minutes).
    pub compression_artifact_remaining_min: f64,
    /// Time of last recorded sample.
    pub last_sample_time_min: f64,
    /// Last published sensor reading.
    pub last_reading_mg_dl: f64,
}

/// CGM simulation engine.
pub struct CgmSensor {
    pub params: CgmParams,
    pub state: CgmState,
    rng: StdRng,
}

impl CgmSensor {
    /// Create a new CGM sensor with optional random seed for reproducibility.
    pub fn new(params: CgmParams, initial_glucose: f64, seed: Option<u64>) -> Self {
        let rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        let state = CgmState {
            interstitial_glucose_mg_dl: initial_glucose,
            ar1_noise_val: 0.0,
            compression_artifact_remaining_min: 0.0,
            last_sample_time_min: -999.0,
            last_reading_mg_dl: initial_glucose,
        };

        Self { params, state, rng }
    }

    /// Advance interstitial diffusion dynamics over step `dt` (min) and optionally emit a new sample.
    ///
    /// # Arguments
    /// * `t`: Current simulation time (minutes).
    /// * `dt`: Integration time step (minutes).
    /// * `plasma_glucose_mg_dl`: Ground truth reference plasma glucose concentration.
    ///
    /// Returns `Some(reading_mg_dl)` if a discrete sensor sample was captured at this time, else `None`.
    pub fn step(&mut self, t: f64, dt: f64, plasma_glucose_mg_dl: f64) -> Option<f64> {
        let p = &self.params;
        let s = &mut self.state;

        // 1. Interstitial glucose diffusion ODE:
        // dG_sub / dt = (G_plasma - G_sub) / tau_cgm
        let dg_sub =
            ((plasma_glucose_mg_dl - s.interstitial_glucose_mg_dl) / p.tau_interstitial_min) * dt;
        s.interstitial_glucose_mg_dl = (s.interstitial_glucose_mg_dl + dg_sub).max(10.0);

        // 2. Check if a periodic discrete sample occurs
        let time_since_last_sample = t - s.last_sample_time_min;
        if time_since_last_sample < p.sampling_interval_min - 1e-6 {
            return None;
        }

        s.last_sample_time_min = t;

        // 3. Update AR(1) colored electrochemical sensor noise
        // epsilon_k = phi * epsilon_{k-1} + w_k, where w_k ~ N(0, sigma^2)
        let normal = Normal::new(0.0, p.noise_std_mg_dl).unwrap();
        let white_noise = normal.sample(&mut self.rng);
        s.ar1_noise_val = p.ar1_phi * s.ar1_noise_val + white_noise;

        // 4. Linear calibration drift over time
        let days_elapsed = t / 1440.0;
        let drift_multiplier = 1.0 + (p.drift_percent_per_day / 100.0) * days_elapsed;

        // 5. Stochastic compression low artifact (mechanical pressure on sensor tissue reduces local perfusion)
        if s.compression_artifact_remaining_min > 0.0 {
            s.compression_artifact_remaining_min -= p.sampling_interval_min;
        } else {
            // Chance of triggering a compression low
            let step_prob = (p.compression_low_prob_per_hr / 60.0) * p.sampling_interval_min;
            if self.rng.gen::<f64>() < step_prob {
                s.compression_artifact_remaining_min = self.rng.gen_range(20.0..60.0);
            }
        }

        let compression_drop = if s.compression_artifact_remaining_min > 0.0 {
            // Drop reading by 25 to 50 mg/dL during compression
            35.0
        } else {
            0.0
        };

        // Synthesize raw measured sensor value
        let raw_cgm =
            (s.interstitial_glucose_mg_dl * drift_multiplier) + s.ar1_noise_val - compression_drop;
        let clamped_reading = raw_cgm.clamp(p.min_reading_mg_dl, p.max_reading_mg_dl);

        // Quantize to 1 mg/dL display resolution (standard in commercial CGMs like Dexcom/Freestyle)
        let quantized_reading = clamped_reading.round();
        s.last_reading_mg_dl = quantized_reading;

        Some(quantized_reading)
    }

    /// Read the latest available interstitial glucose concentration (smooth).
    pub fn interstitial_glucose(&self) -> f64 {
        self.state.interstitial_glucose_mg_dl
    }
}
