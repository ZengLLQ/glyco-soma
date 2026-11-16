//! Virtual patient cohort generation for in-silico clinical trials and population simulations.

use crate::types::{PatientCohort, PatientProfile};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, LogNormal, Normal};

/// Virtual patient cohort generator.
pub struct VirtualCohortGenerator {
    rng: StdRng,
}

impl VirtualCohortGenerator {
    /// Construct generator with optional seed.
    pub fn new(seed: Option<u64>) -> Self {
        let rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };
        Self { rng }
    }

    /// Generate an in-silico population cohort with realistic inter-patient variability.
    pub fn generate_cohort(&mut self, size: usize) -> Vec<PatientProfile> {
        let mut cohort = Vec::with_capacity(size);

        for i in 1..=size {
            // Partition cohort into adult (60%), adolescent (25%), pediatric (15%)
            let roll: f64 = self.rng.gen();
            let cohort_type = if roll < 0.60 {
                PatientCohort::Adult
            } else if roll < 0.85 {
                PatientCohort::Adolescent
            } else {
                PatientCohort::Pediatric
            };

            let patient = self.generate_patient(format!("VP-{:03}", i), cohort_type);
            cohort.push(patient);
        }

        cohort
    }

    /// Generate a single virtual patient belonging to a designated cohort.
    pub fn generate_patient(&mut self, id: String, cohort: PatientCohort) -> PatientProfile {
        match cohort {
            PatientCohort::Adult => {
                let weight_dist = Normal::new(72.0, 12.0).unwrap();
                let raw_weight: f64 = weight_dist.sample(&mut self.rng);
                let weight = raw_weight.clamp(45.0, 120.0);

                // Basal insulin scaled with body weight (~0.012 U/kg/hr)
                let basal_mean = 0.012 * weight;
                let basal_dist = Normal::new(basal_mean, basal_mean * 0.18).unwrap();
                let raw_basal: f64 = basal_dist.sample(&mut self.rng);
                let basal = raw_basal.clamp(0.3, 2.5);

                // CIR log-normal distribution (~10 g/U)
                let cir_dist = LogNormal::new(2.30, 0.28).unwrap();
                let raw_cir: f64 = cir_dist.sample(&mut self.rng);
                let cir = raw_cir.clamp(5.0, 25.0);

                // ISF inversely related to weight and basal (rule of 1800: ISF ~ 1800 / TDD)
                let tdd = basal * 24.0 * 2.0; // TDD ~ 2 * basal
                let isf_base = 1800.0 / tdd.max(15.0);
                let isf_dist = Normal::new(isf_base, isf_base * 0.20).unwrap();
                let raw_isf: f64 = isf_dist.sample(&mut self.rng);
                let isf = raw_isf.clamp(20.0, 90.0);

                let tau_dist = Normal::new(55.0, 8.0).unwrap();
                let raw_tau: f64 = tau_dist.sample(&mut self.rng);
                let tau = raw_tau.clamp(35.0, 80.0);

                PatientProfile {
                    id,
                    cohort,
                    body_weight_kg: weight,
                    basal_rate_u_per_hr: basal,
                    carb_ratio_g_per_u: cir,
                    isf_mg_dl_per_u: isf,
                    target_glucose_mg_dl: 105.0,
                    dia_minutes: 300.0,
                    egp_0: 0.85 * (weight / 70.0),
                    f01_c: 0.80 * (weight / 70.0),
                    renal_threshold_mmol_l: 9.0,
                    vi_l_per_kg: 0.12,
                    vg_l_per_kg: 0.16,
                    tau_insulin_min: tau,
                    dawn_phenomenon_amp: self.rng.gen_range(0.10..0.30),
                }
            }
            PatientCohort::Adolescent => {
                let weight_dist = Normal::new(56.0, 9.0).unwrap();
                let raw_weight: f64 = weight_dist.sample(&mut self.rng);
                let weight = raw_weight.clamp(38.0, 85.0);

                let basal_mean = 0.018 * weight; // adolescent pubertal resistance
                let basal_dist = Normal::new(basal_mean, basal_mean * 0.22).unwrap();
                let raw_basal: f64 = basal_dist.sample(&mut self.rng);
                let basal = raw_basal.clamp(0.5, 3.0);

                let cir_dist = LogNormal::new(2.05, 0.25).unwrap();
                let raw_cir: f64 = cir_dist.sample(&mut self.rng);
                let cir = raw_cir.clamp(4.0, 16.0);

                let tdd = basal * 24.0 * 2.0;
                let isf_base = 1700.0 / tdd.max(20.0);
                let isf_dist = Normal::new(isf_base, isf_base * 0.22).unwrap();
                let raw_isf: f64 = isf_dist.sample(&mut self.rng);
                let isf = raw_isf.clamp(18.0, 60.0);

                PatientProfile {
                    id,
                    cohort,
                    body_weight_kg: weight,
                    basal_rate_u_per_hr: basal,
                    carb_ratio_g_per_u: cir,
                    isf_mg_dl_per_u: isf,
                    target_glucose_mg_dl: 110.0,
                    dia_minutes: 270.0,
                    egp_0: 0.95 * (weight / 55.0),
                    f01_c: 0.90 * (weight / 55.0),
                    renal_threshold_mmol_l: 9.0,
                    vi_l_per_kg: 0.13,
                    vg_l_per_kg: 0.17,
                    tau_insulin_min: self.rng.gen_range(40.0..60.0),
                    dawn_phenomenon_amp: self.rng.gen_range(0.20..0.45),
                }
            }
            PatientCohort::Pediatric => {
                let weight_dist = Normal::new(28.0, 6.0).unwrap();
                let raw_weight: f64 = weight_dist.sample(&mut self.rng);
                let weight = raw_weight.clamp(15.0, 42.0);

                let basal_mean = 0.011 * weight;
                let basal_dist = Normal::new(basal_mean, basal_mean * 0.20).unwrap();
                let raw_basal: f64 = basal_dist.sample(&mut self.rng);
                let basal = raw_basal.clamp(0.15, 0.8);

                let cir_dist = LogNormal::new(2.80, 0.30).unwrap();
                let raw_cir: f64 = cir_dist.sample(&mut self.rng);
                let cir = raw_cir.clamp(10.0, 35.0);

                let tdd = basal * 24.0 * 2.0;
                let isf_base = 1900.0 / tdd.max(8.0);
                let isf_dist = Normal::new(isf_base, isf_base * 0.25).unwrap();
                let raw_isf: f64 = isf_dist.sample(&mut self.rng);
                let isf = raw_isf.clamp(45.0, 140.0);

                PatientProfile {
                    id,
                    cohort,
                    body_weight_kg: weight,
                    basal_rate_u_per_hr: basal,
                    carb_ratio_g_per_u: cir,
                    isf_mg_dl_per_u: isf,
                    target_glucose_mg_dl: 110.0,
                    dia_minutes: 240.0,
                    egp_0: 0.50 * (weight / 28.0),
                    f01_c: 0.55 * (weight / 28.0),
                    renal_threshold_mmol_l: 9.5,
                    vi_l_per_kg: 0.14,
                    vg_l_per_kg: 0.18,
                    tau_insulin_min: self.rng.gen_range(32.0..50.0),
                    dawn_phenomenon_amp: self.rng.gen_range(0.08..0.22),
                }
            }
        }
    }
}
