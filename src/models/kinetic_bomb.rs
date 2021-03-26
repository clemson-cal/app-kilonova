use std::f64::consts::PI;
use serde::{Serialize, Deserialize};
use crate::traits::InitialModel;
use crate::physics::{AnyPrimitive, LIGHT_SPEED};

/**
 * Explosion in a horizontally stratified external medium
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct KineticBomb {
    pub external_medium_density: f64,
    pub launch_radius: f64,
    pub shell_thickness: f64,
    pub kinetic_energy: f64,
    pub shell_mass: f64,
}

// ============================================================================
impl InitialModel for KineticBomb {

    fn validate(&self) -> anyhow::Result<()> {
        Ok(())
    }

    fn primitive_at(&self, coordinate: (f64, f64), time: f64) -> AnyPrimitive {
        
        let (r, _q) = coordinate;

        let v = (2.0 * self.kinetic_energy / self.shell_mass).powf(1.0/2.0) * LIGHT_SPEED;

        if r < v * (time - (self.shell_thickness / v)) || r > v * time {
            let d0 = self.external_medium_density;
            let d = d0 * (r / self.launch_radius).powi(-2);
            let p = d * 1e-3;
        
            AnyPrimitive{
                velocity_r: 0.0,
                velocity_q: 0.0,
                mass_density: d,
                gas_pressure: p,
            }
        } else if r > v * (time - (self.shell_thickness / v)) && r < v * time {
            let vol = 4.0 / 3.0 * PI * ((v * time).powi(3) - (v * (time - (self.shell_thickness / v))).powi(3));
            let d = self.shell_mass / vol;
            let p = d * 1e-3;

            AnyPrimitive{
                velocity_r: v,
                velocity_q: 0.0,
                mass_density: d,
                gas_pressure: p,
            }
        } else {
            AnyPrimitive{
                velocity_r: 0.0,
                velocity_q: 0.0,
                mass_density: 1.0,
                gas_pressure: 0.001,
            }

        }
        
    }

    fn scalar_at(&self, _coordinate: (f64, f64), _time: f64) -> f64 {
        0.0
    }
}
