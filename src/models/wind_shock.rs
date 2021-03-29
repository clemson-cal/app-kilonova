use crate::physics::{AnyPrimitive, LIGHT_SPEED};
use crate::traits::InitialModel;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

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

    /// Rate of outflow of the flare
    pub flare_outflow_rate: f64,

    /// Four velocity of flare
    pub flare_gamma_beta: f64,

    /// Four velocity of wind
    pub wind_gamma_beta: f64,

    /// Four velocity of wind after shock
    pub post_shock_gamma_beta: f64,

    /// Wind pressure
    pub wind_pressure: f64,

    /// Pressure after shock
    pub post_shock_pressure: f64,

    /// Shock location coordinate
    pub shock_location: f64,

    /// Flare time
    pub flare_time: f64,
}

// ============================================================================
impl InitialModel for WindShock {
    fn validate(&self) -> anyhow::Result<()> {
        if self.wind_gamma_beta < 0.0 {
            anyhow::bail!("the wind four-velocity must be positive")
        }
        Ok(())
    }

    fn primitive_at(&self, coordinate: (f64, f64), t: f64) -> AnyPrimitive {
        // u: gamma-beta-c
        // v: beta-c
        // rho: comoving rest-mass density
        // Mdot = 4 pi r^2 rho u c

        if (t >= self.flare_time) & (t < (self.flare_time + 0.1)) {
            let r = coordinate.0;
            let u = self.flare_gamma_beta;
            let rho = self.flare_outflow_rate / (4.0 * PI * r * r * u * LIGHT_SPEED);
            //let rho = (1.0 - n) / (0.1) * (t - self.flare_time - 0.1) + 1.0;
            let p = rho * UNIFORM_TEMPERATURE;
            AnyPrimitive {
                velocity_r: u,
                velocity_q: 0.0,
                mass_density: rho,
                gas_pressure: p,
            }
        } else {
            let r = coordinate.0;
            let u = if r < self.shock_location {
                self.wind_gamma_beta
            } else {
                self.post_shock_gamma_beta
            };
            let rho = self.wind_mass_outflow_rate / (4.0 * PI * r * r * u * LIGHT_SPEED);
            let p = if r < self.shock_location {
                self.wind_pressure
            } else {
                self.post_shock_pressure
            };
            AnyPrimitive {
                velocity_r: u,
                velocity_q: 0.0,
                mass_density: rho,
                gas_pressure: p,
            }
        }
    }

    fn scalar_at(&self, _coordinate: (f64, f64), _t: f64) -> f64 {
        0.0
    }
}
