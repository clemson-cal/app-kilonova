use crate::physics::{AnyPrimitive, LIGHT_SPEED};
use crate::traits::InitialModel;
use serde::{Deserialize, Serialize};

static UNIFORM_TEMPERATURE: f64 = 1e-6;

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

    /// Shock location coordinate
    pub shock_location: f64,
}

// ============================================================================
impl InitialModel for WindShock {
    fn validate(&self) -> anyhow::Result<()> {
        if self.wind_gamma_beta < 0.0 {
            anyhow::bail!("the wind four-velocity must be positive")
        }
        Ok(())
    }

    fn primitive_at(&self, coordinate: (f64, f64), _t: f64) -> AnyPrimitive {
        // u: gamma-beta-c
        // v: beta-c
        // rho: comoving rest-mass density
        // Mdot = 4 pi r^2 rho u c

        let c = 1.0;
        let r = coordinate.0;
        let u = if r < self.shock_location {
            self.wind_gamma_beta * c
        } else {
            self.post_shock_gamma_beta * c
        };
        let rho = self.wind_mass_outflow_rate / (4.0 * PI * r * r * u);

        AnyPrimitive {
            velocity_r: u,
            velocity_q: 0.0,
            mass_density: rho,
            gas_pressure: rho * UNIFORM_TEMPERATURE,
        }
    }

    fn scalar_at(&self, _coordinate: (f64, f64), _t: f64) -> f64 {
        0.0
    }
}
