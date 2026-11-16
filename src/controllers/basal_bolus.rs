//! Standard clinical open-loop basal-bolus therapy with Insulin On Board (IOB) decay tracking.

use crate::types::{HormoneInfusion, PatientProfile};

/// Delivered bolus record for Insulin On Board (IOB) tracking.
#[derive(Debug, Clone, Copy)]
pub struct BolusRecord {
    pub time_min: f64,
    pub amount_u: f64,
}

/// Traditional open-loop insulin pump and multiple daily injections (MDI) calculator.
#[derive(Debug, Clone)]
pub struct BasalBolusController {
    pub patient: PatientProfile,
    pub bolus_history: Vec<BolusRecord>,
}

impl BasalBolusController {
    /// Create new basal-bolus controller for a given patient profile.
    pub fn new(patient: PatientProfile) -> Self {
        Self {
            patient,
            bolus_history: Vec::new(),
        }
    }

    /// Calculate active Insulin On Board ($IOB$) at time `current_time_min` using bilinear/exponential decay.
    pub fn calculate_iob(&self, current_time_min: f64) -> f64 {
        let dia = self.patient.dia_minutes;
        let mut total_iob = 0.0;

        for record in &self.bolus_history {
            let elapsed = current_time_min - record.time_min;
            if elapsed >= 0.0 && elapsed < dia {
                // Standard bilinear decay curve (Mudaliar / Scheiner clinical formula)
                // Fraction remaining: 1.0 - (elapsed / dia) * (1.0 + 0.2 * (1.0 - elapsed / dia))
                let tau = elapsed / dia;
                let fraction_remaining = (1.0 - tau).powi(2);
                total_iob += record.amount_u * fraction_remaining;
            }
        }

        total_iob
    }

    /// Calculate recommended prandial meal bolus:
    ///
    /// $$ B_{meal} = \frac{\text{Carbohydrates (g)}}{CIR} $$
    pub fn calculate_meal_bolus(&self, carbs_grams: f64) -> f64 {
        if self.patient.carb_ratio_g_per_u > 0.0 {
            carbs_grams / self.patient.carb_ratio_g_per_u
        } else {
            0.0
        }
    }

    /// Calculate recommended correction bolus accounting for Insulin On Board ($IOB$):
    ///
    /// $$ B_{corr} = \max\left(0, \frac{G_{meas} - G_{target}}{ISF} - IOB\right) $$
    pub fn calculate_correction_bolus(
        &self,
        current_time_min: f64,
        current_glucose_mg_dl: f64,
    ) -> f64 {
        let target = self.patient.target_glucose_mg_dl;
        if current_glucose_mg_dl <= target {
            return 0.0;
        }

        let isf = self.patient.isf_mg_dl_per_u;
        let raw_correction = (current_glucose_mg_dl - target) / isf;
        let iob = self.calculate_iob(current_time_min);

        (raw_correction - iob).max(0.0)
    }

    /// Register delivery of a meal or correction bolus.
    pub fn deliver_bolus(&mut self, time_min: f64, amount_u: f64) {
        if amount_u > 0.0 {
            self.bolus_history.push(BolusRecord { time_min, amount_u });
        }
    }

    /// Step the basal controller for one time step, delivering nominal basal rate.
    pub fn step(&self) -> HormoneInfusion {
        HormoneInfusion {
            insulin_basal_u_per_hr: self.patient.basal_rate_u_per_hr,
            insulin_bolus_u: 0.0,
            glucagon_basal_ug_per_hr: 0.0,
            glucagon_bolus_ug: 0.0,
        }
    }
}
