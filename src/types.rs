//! Fundamental physiological data types, physical units, and state structures for GlycoSoma.

use serde::{Deserialize, Serialize};

/// Conversion factor: 1 mmol/L of glucose = 18.0182 mg/dL.
pub const MMOL_L_TO_MG_DL: f64 = 18.0182;

/// Conversion factor: 1 mg/dL of glucose = 1 / 18.0182 mmol/L.
pub const MG_DL_TO_MMOL_L: f64 = 1.0 / 18.0182;

/// Conversion factor: 1 Unit of rapid-acting insulin = 6000 pmol.
pub const INSULIN_U_TO_PMOL: f64 = 6000.0;

/// Conversion factor: 1 pmol of insulin = 1 / 6000 Units.
pub const INSULIN_PMOL_TO_U: f64 = 1.0 / 6000.0;

/// Standard normal plasma volume factor (L/kg).
pub const DEFAULT_PLASMA_VOLUME_FACTOR: f64 = 0.045;

/// Clinical glucose units.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GlucoseUnit {
    /// Milligrams per deciliter (standard in US, Japan, France).
    MgDl,
    /// Millimoles per liter (standard in UK, Canada, Australia, Europe).
    MmolL,
}

impl GlucoseUnit {
    /// Convert a value from this unit to mg/dL.
    pub fn to_mg_dl(self, val: f64) -> f64 {
        match self {
            GlucoseUnit::MgDl => val,
            GlucoseUnit::MmolL => val * MMOL_L_TO_MG_DL,
        }
    }

    /// Convert a value from mg/dL to this unit.
    pub fn from_mg_dl(self, val_mg_dl: f64) -> f64 {
        match self {
            GlucoseUnit::MgDl => val_mg_dl,
            GlucoseUnit::MmolL => val_mg_dl * MG_DL_TO_MMOL_L,
        }
    }
}

/// Simulation execution settings.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationConfig {
    /// Baseline discrete integration step size in minutes (e.g. 1.0 min or 0.1 min).
    pub dt: f64,
    /// Total simulation duration in minutes (e.g. 1440.0 for 24 hours).
    pub total_time: f64,
    /// Adaptive RKF45 local truncation error tolerance.
    pub rkf_tolerance: f64,
    /// Minimum allowed time step for adaptive integrators (minutes).
    pub min_dt: f64,
    /// Maximum allowed time step for adaptive integrators (minutes).
    pub max_dt: f64,
    /// Continuous Glucose Monitor (CGM) sampling period in minutes (typically 5.0 min).
    pub cgm_sampling_interval: f64,
    /// Random seed for stochastic noise repeatability.
    pub random_seed: Option<u64>,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            dt: 1.0,
            total_time: 1440.0, // 24 hours
            rkf_tolerance: 1e-5,
            min_dt: 1e-4,
            max_dt: 5.0,
            cgm_sampling_interval: 5.0,
            random_seed: Some(42),
        }
    }
}

/// Virtual patient anthropometric and metabolic characteristics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PatientProfile {
    /// Patient identifier or name.
    pub id: String,
    /// Patient age cohort.
    pub cohort: PatientCohort,
    /// Body weight in kilograms (kg).
    pub body_weight_kg: f64,
    /// Baseline basal insulin infusion rate (Units/hour).
    pub basal_rate_u_per_hr: f64,
    /// Carbohydrate-to-Insulin Ratio (CIR, grams of carbohydrate per Unit of insulin).
    pub carb_ratio_g_per_u: f64,
    /// Insulin Sensitivity Factor (ISF, drop in glucose in mg/dL per Unit of insulin).
    pub isf_mg_dl_per_u: f64,
    /// Target glucose setpoint in mg/dL (typically 100.0 to 110.0 mg/dL).
    pub target_glucose_mg_dl: f64,
    /// Duration of Insulin Action (DIA) in minutes (typically 240 to 360 min).
    pub dia_minutes: f64,
    /// Endogenous glucose production baseline at zero remote insulin (mmol/min).
    pub egp_0: f64,
    /// Non-insulin-dependent glucose uptake (CNS/erythrocytes) (mmol/min).
    pub f01_c: f64,
    /// Renal glucose clearance threshold in mmol/L (typically 9.0 mmol/L = ~162 mg/dL).
    pub renal_threshold_mmol_l: f64,
    /// Insulin distribution volume (L/kg).
    pub vi_l_per_kg: f64,
    /// Glucose distribution volume (L/kg).
    pub vg_l_per_kg: f64,
    /// Subcutaneous insulin absorption time constant (minutes).
    pub tau_insulin_min: f64,
    /// Circadian dawn phenomenon amplitude (fractional increase in hepatic output 04:00-08:00).
    pub dawn_phenomenon_amp: f64,
}

/// Age categories with distinct metabolic parameter distributions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PatientCohort {
    /// Adult patient (18–65 years).
    Adult,
    /// Adolescent patient (12–18 years), exhibiting pubertal insulin resistance and hormonal surges.
    Adolescent,
    /// Pediatric patient (4–11 years), high insulin sensitivity and rapid carbohydrate absorption.
    Pediatric,
}

impl Default for PatientProfile {
    fn default() -> Self {
        Self::standard_adult()
    }
}

impl PatientProfile {
    /// Create a representative standard adult Type 1 Diabetes patient.
    pub fn standard_adult() -> Self {
        Self {
            id: "Adult-001".to_string(),
            cohort: PatientCohort::Adult,
            body_weight_kg: 70.0,
            basal_rate_u_per_hr: 0.85,
            carb_ratio_g_per_u: 10.0,
            isf_mg_dl_per_u: 45.0,
            target_glucose_mg_dl: 105.0,
            dia_minutes: 300.0, // 5 hours
            egp_0: 1.05,        // mmol/min
            f01_c: 0.45,        // mmol/min
            renal_threshold_mmol_l: 9.0,
            vi_l_per_kg: 0.12,
            vg_l_per_kg: 0.16,
            tau_insulin_min: 55.0,
            dawn_phenomenon_amp: 0.20,
        }
    }

    /// Create an adolescent patient with increased hormonal resistance.
    pub fn adolescent() -> Self {
        Self {
            id: "Adolescent-001".to_string(),
            cohort: PatientCohort::Adolescent,
            body_weight_kg: 55.0,
            basal_rate_u_per_hr: 1.10,
            carb_ratio_g_per_u: 7.5,
            isf_mg_dl_per_u: 30.0,
            target_glucose_mg_dl: 110.0,
            dia_minutes: 270.0,
            egp_0: 1.15,
            f01_c: 0.50,
            renal_threshold_mmol_l: 9.0,
            vi_l_per_kg: 0.13,
            vg_l_per_kg: 0.17,
            tau_insulin_min: 48.0,
            dawn_phenomenon_amp: 0.35,
        }
    }

    /// Create a pediatric patient with high insulin sensitivity.
    pub fn pediatric() -> Self {
        Self {
            id: "Pediatric-001".to_string(),
            cohort: PatientCohort::Pediatric,
            body_weight_kg: 28.0,
            basal_rate_u_per_hr: 0.35,
            carb_ratio_g_per_u: 18.0,
            isf_mg_dl_per_u: 85.0,
            target_glucose_mg_dl: 110.0,
            dia_minutes: 240.0,
            egp_0: 0.65,
            f01_c: 0.28,
            renal_threshold_mmol_l: 9.5,
            vi_l_per_kg: 0.14,
            vg_l_per_kg: 0.18,
            tau_insulin_min: 40.0,
            dawn_phenomenon_amp: 0.15,
        }
    }

    /// Compute total glucose distribution volume in liters ($V_G = V_{G,factor} \times W$).
    pub fn total_glucose_volume_l(&self) -> f64 {
        self.vg_l_per_kg * self.body_weight_kg
    }

    /// Compute total insulin distribution volume in liters ($V_I = V_{I,factor} \times W$).
    pub fn total_insulin_volume_l(&self) -> f64 {
        self.vi_l_per_kg * self.body_weight_kg
    }
}

/// Scheduled or unannounced carbohydrate ingestion event.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MealEvent {
    /// Ingestion start time in simulation minutes.
    pub start_time_min: f64,
    /// Total grams of carbohydrates (CHO, grams).
    pub carbs_grams: f64,
    /// Meal duration in minutes (e.g. 15.0 to 30.0 min).
    pub duration_min: f64,
    /// Glycemic Index category or fractional absorption speed (0.3 = low GI, 1.0 = rapid simple sugar).
    pub glycemic_index: f64,
    /// Optional fat and protein content in grams (causing delayed gastric emptying).
    pub fat_grams: f64,
    pub protein_grams: f64,
}

impl MealEvent {
    /// Create a standard rapid meal.
    pub fn standard(start_time_min: f64, carbs_grams: f64) -> Self {
        Self {
            start_time_min,
            carbs_grams,
            duration_min: 20.0,
            glycemic_index: 0.75,
            fat_grams: 10.0,
            protein_grams: 15.0,
        }
    }

    /// Create a high-fat slow meal (e.g. pizza).
    pub fn slow_high_fat(start_time_min: f64, carbs_grams: f64) -> Self {
        Self {
            start_time_min,
            carbs_grams,
            duration_min: 45.0,
            glycemic_index: 0.50,
            fat_grams: 35.0,
            protein_grams: 30.0,
        }
    }
}

/// Physical exercise exertion event modifying metabolic consumption.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExerciseEvent {
    /// Exercise onset time in simulation minutes.
    pub start_time_min: f64,
    /// Exercise duration in minutes.
    pub duration_min: f64,
    /// Metabolic Equivalent of Task (METs: 1.0 = resting, 4.0 = moderate walk, 8.0 = intense running).
    pub intensity_mets: f64,
    /// Fraction of heart rate reserve (0.0 to 1.0).
    pub hr_reserve_fraction: f64,
}

impl ExerciseEvent {
    /// Moderate aerobic activity (e.g. brisk walking or cycling at 5 METs).
    pub fn moderate_aerobic(start_time_min: f64, duration_min: f64) -> Self {
        Self {
            start_time_min,
            duration_min,
            intensity_mets: 5.0,
            hr_reserve_fraction: 0.60,
        }
    }

    /// Vigorous high-intensity interval activity (8 METs).
    pub fn vigorous(start_time_min: f64, duration_min: f64) -> Self {
        Self {
            start_time_min,
            duration_min,
            intensity_mets: 8.0,
            hr_reserve_fraction: 0.85,
        }
    }
}

/// Dual-hormone actuator delivery commands.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct HormoneInfusion {
    /// Continuous basal insulin infusion rate (Units/hour).
    pub insulin_basal_u_per_hr: f64,
    /// Discrete bolus insulin delivery (Units).
    pub insulin_bolus_u: f64,
    /// Continuous glucagon infusion rate (micrograms/hour).
    pub glucagon_basal_ug_per_hr: f64,
    /// Discrete glucagon rescue micro-bolus (micrograms).
    pub glucagon_bolus_ug: f64,
}

impl HormoneInfusion {
    /// Convert total insulin delivery rate over time step `dt` (min) into Units.
    pub fn total_insulin_u_in_step(&self, dt_min: f64) -> f64 {
        (self.insulin_basal_u_per_hr / 60.0) * dt_min + self.insulin_bolus_u
    }

    /// Convert total glucagon delivery rate over time step `dt` (min) into micrograms.
    pub fn total_glucagon_ug_in_step(&self, dt_min: f64) -> f64 {
        (self.glucagon_basal_ug_per_hr / 60.0) * dt_min + self.glucagon_bolus_ug
    }
}

/// Complete multi-compartmental physiological state vector at a specific time point.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetabolicState {
    /// Elapsed simulation time in minutes.
    pub time_min: f64,
    /// Plasma glucose concentration in mg/dL.
    pub glucose_plasma_mg_dl: f64,
    /// Subcutaneous interstitial glucose concentration in mg/dL (sensed by CGM).
    pub glucose_interstitial_mg_dl: f64,
    /// Plasma insulin concentration in mU/L.
    pub insulin_plasma_mu_l: f64,
    /// First subcutaneous insulin depot compartment ($S_1$, pmol).
    pub s1_subcut_insulin: f64,
    /// Second subcutaneous insulin depot compartment ($S_2$, pmol).
    pub s2_subcut_insulin: f64,
    /// Remote insulin action on glucose transport ($x_1$, 1/min).
    pub x1_remote_transport: f64,
    /// Remote insulin action on glucose disposal ($x_2$, 1/min).
    pub x2_remote_disposal: f64,
    /// Remote insulin action on endogenous glucose production ($x_3$, unitless).
    pub x3_remote_egp: f64,
    /// Accessible glucose mass ($Q_1$, mmol).
    pub q1_accessible_glucose_mmol: f64,
    /// Non-accessible glucose mass ($Q_2$, mmol).
    pub q2_nonaccessible_glucose_mmol: f64,
    /// Solid phase stomach carbohydrate contents (grams).
    pub stomach_solid_g: f64,
    /// Liquid phase stomach carbohydrate contents (grams).
    pub stomach_liquid_g: f64,
    /// Intestinal carbohydrate absorption compartment (grams).
    pub intestine_g: f64,
    /// Plasma glucagon concentration in pg/mL.
    pub glucagon_plasma_pg_ml: f64,
    /// Liver glycogen stores (fraction of capacity 0.0 to 1.0).
    pub liver_glycogen_fraction: f64,
    /// Active muscle non-insulin-mediated glucose uptake rate (mmol/min).
    pub active_exercise_uptake_mmol_min: f64,
}

impl MetabolicState {
    /// Initialize homeostatic basal state for a given patient profile.
    pub fn initial_basal(profile: &PatientProfile) -> Self {
        let vg = profile.total_glucose_volume_l();
        let target_mmol_l = profile.target_glucose_mg_dl * MG_DL_TO_MMOL_L;
        let q1_init = target_mmol_l * vg;
        let q2_init = q1_init * 0.824;

        // Steady state insulin calculations
        let basal_pmol_min = (profile.basal_rate_u_per_hr * INSULIN_U_TO_PMOL) / 60.0;
        let tau_s = profile.tau_insulin_min;
        let s1_init = basal_pmol_min * tau_s;
        let s2_init = s1_init;
        let vi = profile.total_insulin_volume_l();
        let i_init_pmol_l = (basal_pmol_min / vi) / 0.138;
        let i_init_mu_l = i_init_pmol_l / 6.0;

        let x1_init = (5.12e-6 / 0.006) * i_init_pmol_l;
        let x2_init = (8.22e-6 / 0.060) * i_init_pmol_l;
        let x3_init = (2.60e-6 / 0.030) * i_init_pmol_l;

        Self {
            time_min: 0.0,
            glucose_plasma_mg_dl: profile.target_glucose_mg_dl,
            glucose_interstitial_mg_dl: profile.target_glucose_mg_dl,
            insulin_plasma_mu_l: i_init_mu_l,
            s1_subcut_insulin: s1_init,
            s2_subcut_insulin: s2_init,
            x1_remote_transport: x1_init,
            x2_remote_disposal: x2_init,
            x3_remote_egp: x3_init,
            q1_accessible_glucose_mmol: q1_init,
            q2_nonaccessible_glucose_mmol: q2_init,
            stomach_solid_g: 0.0,
            stomach_liquid_g: 0.0,
            intestine_g: 0.0,
            glucagon_plasma_pg_ml: 65.0,
            liver_glycogen_fraction: 0.90,
            active_exercise_uptake_mmol_min: 0.0,
        }
    }
}
