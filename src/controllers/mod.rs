//! Automated artificial pancreas controllers: Model Predictive Control, Dual-Hormone PID, PLGS, and Basal-Bolus therapy.

pub mod basal_bolus;
pub mod mpc;
pub mod pid;
pub mod plgs;

pub use basal_bolus::{BasalBolusController, BolusRecord};
pub use mpc::{ModelPredictiveController, MpcConfig};
pub use pid::{PidConfig, PidController};
pub use plgs::{PlgsConfig, PlgsController, PlgsState};

/// Unified controller modality selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ControllerMode {
    /// Traditional open-loop scheduled basal with manual meal boluses.
    OpenLoopBasalBolus,
    /// Automated Predictive Low Glucose Suspend (PLGS) safety mechanism.
    PredictiveLowGlucoseSuspend,
    /// Dual-hormone feedback PID controller.
    DualHormonePid,
    /// Constrained Model Predictive Control (MPC) with asymmetric cost functions.
    ModelPredictiveControl,
}

impl ControllerMode {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            ControllerMode::OpenLoopBasalBolus => "Open-Loop Basal-Bolus",
            ControllerMode::PredictiveLowGlucoseSuspend => "Predictive Low Glucose Suspend (PLGS)",
            ControllerMode::DualHormonePid => "Dual-Hormone PID",
            ControllerMode::ModelPredictiveControl => "Constrained MPC",
        }
    }
}
