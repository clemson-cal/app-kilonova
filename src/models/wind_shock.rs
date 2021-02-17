use crate::physics::AgnosticPrimitive;
use crate::traits::InitialModel;
use serde::{Deserialize, Serialize};

static UNIFORM_TEMPERATURE: f64 = 1e-3;

/**
 * Jet propagating through a kilonova debris cloud and surrounding relativistic
 * envelop
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WindShock {
    /// Rate of outflow of the wind
    pub wind_mass_outflow_rate: f64,

    /// Four velocity of wind
    pub wind_gamma_beta: f64,

    /// Four velocity of wind after shock
    pub post_shock_gamma_beta: f64,
}

// ============================================================================
impl InitialModel for WindShock {
    fn validate(&self) -> anyhow::Result<()> {
        Ok(())
    }

    fn primitive_at(&self, _coordinate: (f64, f64), _t: f64) -> AgnosticPrimitive {
        AgnosticPrimitive {
            velocity_r: self.wind_gamma_beta,
            velocity_q: 0.0,
            mass_density: 1.0,
            gas_pressure: UNIFORM_TEMPERATURE,
        }
    }

    fn scalar_at(&self, _coordinate: (f64, f64), _t: f64) -> f64 {
        0.0
    }
}
