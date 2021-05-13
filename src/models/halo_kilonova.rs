use std::sync::{Arc, Mutex};
use crate::lookup_table_v2::LookupTable;
use crate::physics::{AnyPrimitive, LIGHT_SPEED};
use crate::traits::InitialModel;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

const UNIFORM_TEMPERATURE: f64 = 1e-3;

/**
 * Explosion in a horizontally stratified external medium
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HaloKilonova {
    pub altitude: f64,
    pub launch_radius: f64,
    pub shell_thickness: f64,
    pub kinetic_energy: f64,
    pub shell_mass: f64,
    pub radial_distance: f64,
    pub initial_data_table: Option<String>,

    #[serde(skip)]
    pub lookup_table: Arc<Mutex<Option<LookupTable<3>>>>
}

// ============================================================================
impl HaloKilonova {
    fn shell_extent(&self, t: f64) -> std::ops::Range<f64> {
        let r_outer_shell_surface = self.launch_radius + self.shell_velocity() * t;
        let r_inner_shell_surface =
            self.launch_radius + self.shell_velocity() * (t - self.shell_duration());
        r_inner_shell_surface..r_outer_shell_surface
    }

    fn shell_velocity(&self) -> f64 {
        (2.0 * self.kinetic_energy / self.shell_mass).sqrt()
    }

    fn shell_duration(&self) -> f64 {
        self.shell_thickness / self.shell_velocity()
    }

    fn require_lookup_table(&self) {
        let mut self_table = self.lookup_table.as_ref().lock().unwrap();

        if self_table.is_none() {
            let filename = self.initial_data_table.as_ref().unwrap();
            let table = LookupTable::<3>::from_ascii_file(&filename).unwrap();
            *self_table = Some(table);
        }
    }
}

// ============================================================================
impl InitialModel for HaloKilonova {
    fn validate(&self) -> anyhow::Result<()> {
        if self.shell_velocity() > 0.25 * LIGHT_SPEED {
            anyhow::bail! {"
            The shell is moving faster (v/c = {}) than 0.25 c, but
            this problem assumes Newtonian expressions for the
            kinetic energy. Consider reducing the kinetic energy or
            increasing the shell mass.", self.shell_velocity() / LIGHT_SPEED}
        } else if let Some(initial_data_table) = &self.initial_data_table {
            LookupTable::<3>::from_ascii_file(initial_data_table)?;
        }
        Ok(())
    }

    fn primitive_at(&self, coordinate: (f64, f64), t: f64) -> AnyPrimitive {
        let (r, q) = coordinate;
        let z = r * q.cos() + self.altitude;

        if self.shell_extent(t).contains(&r) {
            let mdot = self.shell_mass / self.shell_duration();
            let v = self.shell_velocity();
            let d = mdot / (4.0 * PI * r * r * v);
            let p = d * UNIFORM_TEMPERATURE;

            AnyPrimitive {
                velocity_r: v / LIGHT_SPEED,
                velocity_q: 0.0,
                mass_density: d,
                gas_pressure: p,
            }
        } else if self.initial_data_table.is_some() {
            self.require_lookup_table();
            let table_borrow = self.lookup_table.as_ref().lock().unwrap();
            let table = table_borrow.as_ref().unwrap();
            let sample = table.sample(z);
            let d = sample[2];
            let p = sample[1] / LIGHT_SPEED / LIGHT_SPEED;

            AnyPrimitive {
                velocity_r: 0.0,
                velocity_q: 0.0,
                mass_density: d,
                gas_pressure: p,
            }

        } else {
            panic!("the halo kilonova setup requires z > 0.0");
        }
    }

    fn scalar_at(&self, coordinate: (f64, f64), t: f64) -> f64 {
        let (r, _q) = coordinate;
        if self.shell_extent(t).contains(&r) {
            1.0
        } else {
            0.0
        }
    }
}
