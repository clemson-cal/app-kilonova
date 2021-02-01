use std::f64::consts::PI;
use serde::{Serialize, Deserialize};
use crate::traits::InitialModel;
use crate::physics::{AgnosticPrimitive, LIGHT_SPEED};




/**
 * Explosion in a horizontally stratified external medium
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HaloKilonova {
    pub altitude: f64,
    pub external_medium_density: f64,
    pub launch_radius: f64,
    pub ejecta_mdot: f64,
    pub ejecta_gamma_beta: f64,
    pub engine_duration: f64,
}




// ============================================================================
impl InitialModel for HaloKilonova {

    fn validate(&self) -> anyhow::Result<()> {
        Ok(())
    }

    fn primitive_at(&self, coordinate: (f64, f64), time: f64) -> AgnosticPrimitive {
        
        let (r, q) = coordinate;

        if r > self.launch_radius {
            let z = r * q.cos();
            let z0 = -self.altitude;
            let d0 = self.external_medium_density;
            let d = d0 * if z < z0 {
                d0
            } else {
                d0 * f64::exp(-f64::powf((z - z0) / z0, 2.0))
            };
            let p = d * 1e-3;
        
            AgnosticPrimitive{
                velocity_r: 0.0,
                velocity_q: 0.0,
                mass_density: d,
                gas_pressure: p,
            }
        } else if time < self.engine_duration {
            let d = self.ejecta_mdot / 4.0 / PI / r / r / self.ejecta_gamma_beta / LIGHT_SPEED;
            let p = d * 1e-3;

            AgnosticPrimitive{
                velocity_r: self.ejecta_gamma_beta,
                velocity_q: 0.0,
                mass_density: d,
                gas_pressure: p,
            }
        } else {
            AgnosticPrimitive{
                velocity_r: 0.0,
                velocity_q: 0.0,
                mass_density: 1.0,
                gas_pressure: 0.1,
            }
        }
        
    }

    fn scalar_at(&self, _coordinate: (f64, f64), _time: f64) -> f64 {
        0.0
    }
}
