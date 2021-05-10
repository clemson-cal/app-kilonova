use crate::lookup_table_v2::LookupTable;
use crate::physics::{AnyPrimitive, LIGHT_SPEED};
use crate::traits::InitialModel;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use std::cell::RefCell;

static UNIFORM_TEMPERATURE: f64 = 1e-6;

/// Jet propagating through a kilonova debris cloud and surrounding
/// relativistic envelop
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WindShock {
    /// Rate of outflow of the wind
    pub wind_mass_outflow_rate: f64,

    /// Four velocity of wind
    pub wind_gamma_beta: f64,

    /// Wind pressure
    pub wind_pressure: f64,

    /// Pressure after shock
    pub post_shock_pressure: f64,

    /// Shock location coordinate
    pub shock_location: f64,

    /// Four velocity of wind after shock
    pub post_shock_gamma_beta: f64,

    /// Rate of outflow of the flare
    #[serde(default)]
    pub flare_outflow_rate: f64,

    /// Four velocity of flare
    #[serde(default)]
    pub flare_gamma_beta: f64,

    /// Flare time
    #[serde(default)]
    pub flare_time: f64,

    /// Flare duration
    #[serde(default)]
    pub flare_duration: f64,

    /// Initial data table. This field is optional. If it's given a value, it
    /// must be the relative path to an ASCII table of initial data for a
    /// wind. The table columns are expected to be (radius [cm], gamma-beta,
    /// mass density [g / cm^3], specific enthalpy [cm^2 / s^2]). If given,
    /// the above parameters are ignored, except for the ones starting with
    /// `flare`.
    pub initial_data_table: Option<String>,

    #[serde(skip)]
    pub lookup_table: RefCell<Option<LookupTable<4>>>,
}

impl WindShock {
    fn require_lookup_table(&self) {
        if self.lookup_table.borrow().is_none() {
            let filename = self.initial_data_table.as_ref().unwrap();
            let table = LookupTable::<4>::from_ascii_file(&filename).unwrap();
            *self.lookup_table.borrow_mut() = Some(table);
        }
    }
}

// ============================================================================
impl InitialModel for WindShock {
    fn validate(&self) -> anyhow::Result<()> {
        if self.wind_gamma_beta < 0.0 {
            anyhow::bail!("the wind four-velocity must be positive")
        } else if let Some(initial_data_table) = &self.initial_data_table {
            LookupTable::<4>::from_ascii_file(initial_data_table)?;
        }
        Ok(())
    }

    fn primitive_at(&self, coordinate: (f64, f64), t: f64) -> AnyPrimitive {
        // u: gamma-beta-c
        // v: beta-c
        // rho: comoving rest-mass density
        // Mdot = 4 pi r^2 rho u c

        if t >= self.flare_time && t < self.flare_time + self.flare_duration {
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
        } else if t >= self.flare_time + 3.0 && t < self.flare_time + 3.0 + self.flare_duration {
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
        } else if t >= self.flare_time + 6.0 && t < self.flare_time + 6.0 + self.flare_duration {
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
        } else if self.initial_data_table.is_some() {
            self.require_lookup_table();
            let table_borrow = self.lookup_table.borrow();
            let table = table_borrow.as_ref().unwrap();
            let sample = table.sample(coordinate.0);
            let u = sample[1];
            let d = sample[2];
            let h = sample[3];
            let mu = h / LIGHT_SPEED / LIGHT_SPEED - 1.0;
            let e = mu / (4.0 / 3.0);
            let p = d * e * (4.0 / 3.0 - 1.0);
            AnyPrimitive {
                velocity_r: u,
                velocity_q: 0.0,
                mass_density: d,
                gas_pressure: p,
            }
        } else {
            todo!("restore evaluation of wind profile which does not rely on a table")
        }
    }

    fn scalar_at(&self, _coordinate: (f64, f64), _t: f64) -> f64 {
        0.0
    }
}
