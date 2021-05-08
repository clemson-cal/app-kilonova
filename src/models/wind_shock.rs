// use crate::ascii_lookup::LookupTable;
use crate::physics::{AnyPrimitive, LIGHT_SPEED};
use crate::traits::InitialModel;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

static UNIFORM_TEMPERATURE: f64 = 1e-6;

/// Jet propagating through a kilonova debris cloud and surrounding
/// relativistic envelop
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WindShock {
    /// Rate of outflow of the wind
    pub wind_mass_outflow_rate: f64,

    /// Rate of outflow of the flare
    #[serde(default)]
    pub flare_outflow_rate: f64,

    /// Four velocity of flare
    #[serde(default)]
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

    /// Initial data _table
    #[serde(default)]
    pub initial_data_table: Option<String>,

    /// Flare time
    #[serde(default)]
    pub flare_time: f64,

    /// Flare duration
    #[serde(default)]
    pub flare_duration: f64,
}

// thread_local! {
//     static DAT: Vec<[f64; 4]> = LookupTable::<4>::from_ascii();
// }

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

        if (t >= self.flare_time) & (t < (self.flare_time + self.flare_duration)) {
            let r = coordinate.0;
            let u = self.flare_gamma_beta;
            let n = self.flare_outflow_rate / (4.0 * PI * r * r * u * LIGHT_SPEED);
            let rho = n * (self.flare_time + self.flare_duration - t) / self.flare_duration;
            let p = rho * UNIFORM_TEMPERATURE;

            AnyPrimitive {
                velocity_r: u,
                velocity_q: 0.0,
                mass_density: rho,
                gas_pressure: p,
            }
        } else if (t >= self.flare_time + 3.0) & (t < (self.flare_time + 3.0 + self.flare_duration)) {
            let r = coordinate.0;
            let u = self.flare_gamma_beta;
            let n = self.flare_outflow_rate / (4.0 * PI * r * r * u * LIGHT_SPEED);
            let rho = n * (self.flare_time + 3.0 + self.flare_duration - t) / self.flare_duration;
            let p = rho * UNIFORM_TEMPERATURE;

            AnyPrimitive {
                velocity_r: u,
                velocity_q: 0.0,
                mass_density: rho,
                gas_pressure: p,
            }
        } else if (t >= self.flare_time + 6.0) & (t < (self.flare_time + 6.0 + self.flare_duration)) {
            let r = coordinate.0;
            let u = self.flare_gamma_beta;
            let n = self.flare_outflow_rate / (4.0 * PI * r * r * u * LIGHT_SPEED);
            let rho = n * (self.flare_time + 6.0 + self.flare_duration - t) / self.flare_duration;
            let p = rho * UNIFORM_TEMPERATURE;

            AnyPrimitive {
                velocity_r: u,
                velocity_q: 0.0,
                mass_density: rho,
                gas_pressure: p,
            }
        } else {
            todo!()
            // let r = coordinate.0;
            // DAT.with(|f| {
            //     let four_velocity = LookupTable { rows: f.to_vec() }.sample(r)[1];
            //     let mass_density = LookupTable { rows: f.to_vec() }.sample(r)[2];
            //     let sp_enthalpy = LookupTable { rows: f.to_vec() }.sample(r)[3];

            //     let u = four_velocity;
            //     let rho = mass_density;
            //     let h = sp_enthalpy;
            //     let mu = h - LIGHT_SPEED * LIGHT_SPEED;
            //     let e = mu / (4.0 / 3.0);
            //     let p = rho * e * (4.0 / 3.0 - 1.0);
            //     AnyPrimitive {
            //         velocity_r: u,
            //         velocity_q: 0.0,
            //         mass_density: rho,
            //         gas_pressure: p,
            //     }
            // })
        }
    }

    fn scalar_at(&self, _coordinate: (f64, f64), _t: f64) -> f64 {
        0.0
    }
}
