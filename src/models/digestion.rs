//! Gastrointestinal carbohydrate digestion, gastric emptying, and intestinal glucose absorption dynamics.
//!
//! # Mechanisms
//! - Two-compartment stomach (solid vs liquid phase transit).
//! - Negative feedback inhibition of gastric emptying by duodenal chyme osmolarity and glucose presence.
//! - Multi-compartment small-intestinal transit and active SGLT-1 / GLUT-2 absorption.
//! - Incretin effect (GLP-1 release proportional to carbohydrate absorption velocity).

use crate::types::MealEvent;

/// Gastrointestinal physiology parameters.
#[derive(Debug, Clone)]
pub struct DigestionParams {
    /// Baseline gastric emptying rate for liquid phase $k_{empt, liquid}$ (1/min, $t_{1/2} \approx 20$ min).
    pub k_empt_liquid: f64,
    /// Solid phase breakdown/trituration rate constant $k_{trit}$ (1/min, $t_{1/2} \approx 45$ min).
    pub k_trit_solid: f64,
    /// Intestinal absorption rate constant $k_{abs}$ (1/min, $t_{1/2} \approx 15$ min).
    pub k_absorption: f64,
    /// Bioavailability of ingested carbohydrates ($f \approx 0.90$ to 0.95).
    pub bioavailability: f64,
    /// Feedback inhibition constant for gastric emptying from small intestine carbohydrate mass (grams).
    pub k_feedback_mass_g: f64,
    /// Incretin GLP-1 secretion sensitivity per gram of absorbed glucose.
    pub glp1_secretion_per_g: f64,
}

impl Default for DigestionParams {
    fn default() -> Self {
        Self {
            k_empt_liquid: 0.035,
            k_trit_solid: 0.015,
            k_absorption: 0.055,
            bioavailability: 0.92,
            k_feedback_mass_g: 25.0,
            glp1_secretion_per_g: 0.08,
        }
    }
}

/// Dynamic state of the digestive tract.
#[derive(Debug, Clone, Default)]
pub struct DigestionState {
    /// Solid carbohydrate mass in stomach ($Q_{sto1}$, grams).
    pub stomach_solid_g: f64,
    /// Liquid/chyme carbohydrate mass in stomach ready for pyloric transit ($Q_{sto2}$, grams).
    pub stomach_liquid_g: f64,
    /// Proximal small intestine carbohydrate mass ($Q_{gut1}$, grams).
    pub gut_proximal_g: f64,
    /// Distal small intestine carbohydrate mass ($Q_{gut2}$, grams).
    pub gut_distal_g: f64,
    /// Current plasma GLP-1 concentration increment (pmol/L).
    pub glp1_increment_pmol_l: f64,
}

/// Digestion and intestinal carbohydrate absorption simulation engine.
#[derive(Debug, Clone)]
pub struct DigestionSystem {
    pub params: DigestionParams,
    pub state: DigestionState,
}

impl DigestionSystem {
    /// Create new digestion system with specified parameters.
    pub fn new(params: DigestionParams) -> Self {
        Self {
            params,
            state: DigestionState::default(),
        }
    }

    /// Add an ingested meal event into the gastric solid compartment.
    pub fn ingest_meal(&mut self, meal: &MealEvent) {
        // High glycemic index meals partition more directly into liquid/readily emptiable chyme
        let liquid_fraction = (meal.glycemic_index * 0.65).clamp(0.20, 0.85);
        let solid_fraction = 1.0 - liquid_fraction;

        let total_carbs = meal.carbs_grams * self.params.bioavailability;
        self.state.stomach_solid_g += total_carbs * solid_fraction;
        self.state.stomach_liquid_g += total_carbs * liquid_fraction;
    }

    /// Advance digestion state across interval `dt` (min).
    ///
    /// Returns the systemic glucose Rate of Appearance ($R_a$, or $U_G$) in mmol/min.
    /// (Note: 1 gram of glucose = 1000 mg / 180.1558 g/mol = 5.5507 mmol).
    pub fn step(&mut self, dt: f64) -> f64 {
        let p = &self.params;
        let s = &mut self.state;

        // Feedback inhibition on gastric emptying from duodenal glucose presence:
        // Inhibits pyloric opening when duodenum is overloaded with nutrients
        let total_gut = s.gut_proximal_g + s.gut_distal_g;
        let feedback_factor = 1.0 / (1.0 + (total_gut / p.k_feedback_mass_g).powi(2));

        // Trituration: solid -> liquid chyme in stomach
        let trituration_flux = p.k_trit_solid * s.stomach_solid_g * dt;
        let actual_trit = trituration_flux.min(s.stomach_solid_g);
        s.stomach_solid_g -= actual_trit;
        s.stomach_liquid_g += actual_trit;

        // Gastric emptying: liquid stomach -> proximal gut
        let emptying_rate = p.k_empt_liquid * feedback_factor;
        let emptying_flux = emptying_rate * s.stomach_liquid_g * dt;
        let actual_empty = emptying_flux.min(s.stomach_liquid_g);
        s.stomach_liquid_g -= actual_empty;
        s.gut_proximal_g += actual_empty;

        // Intestinal transit: proximal gut -> distal gut
        let transit_flux = p.k_absorption * 0.60 * s.gut_proximal_g * dt;
        let actual_transit = transit_flux.min(s.gut_proximal_g);
        s.gut_proximal_g -= actual_transit;
        s.gut_distal_g += actual_transit;

        // Systemic intestinal absorption: flux across gut mucosa into portal circulation
        let absorption_flux_g =
            (p.k_absorption * s.gut_proximal_g * 0.40 + p.k_absorption * s.gut_distal_g) * dt;
        let actual_absorbed_g = absorption_flux_g.min(s.gut_proximal_g + s.gut_distal_g);

        // Deduct absorbed glucose from distal gut
        let deduct_distal = actual_absorbed_g.min(s.gut_distal_g);
        s.gut_distal_g -= deduct_distal;
        let remaining_to_deduct = actual_absorbed_g - deduct_distal;
        s.gut_proximal_g = (s.gut_proximal_g - remaining_to_deduct).max(0.0);

        // Convert absorbed grams per dt into mmol/min
        // 1 g glucose = 5.5507 mmol
        let grams_per_min = actual_absorbed_g / dt.max(1e-6);
        let mmol_per_min = grams_per_min * 5.55075;

        // Incretin GLP-1 secretion kinetics (elimination half-life ~2 min)
        let glp1_clearance = 0.35; // 1/min
        let glp1_production = grams_per_min * p.glp1_secretion_per_g * 10.0;
        let d_glp1 = (glp1_production - glp1_clearance * s.glp1_increment_pmol_l) * dt;
        s.glp1_increment_pmol_l = (s.glp1_increment_pmol_l + d_glp1).max(0.0);

        mmol_per_min
    }

    /// Total carbohydrate mass currently remaining in the entire gastrointestinal tract (grams).
    pub fn total_gut_carbs_g(&self) -> f64 {
        self.state.stomach_solid_g
            + self.state.stomach_liquid_g
            + self.state.gut_proximal_g
            + self.state.gut_distal_g
    }
}
