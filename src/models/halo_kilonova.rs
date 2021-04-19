use std::f64::consts::PI;
use serde::{Serialize, Deserialize};
use crate::traits::InitialModel;
use crate::physics::{AnyPrimitive, LIGHT_SPEED};
use crate::galmod::GalacticModel;
use float_ord::FloatOrd;

const UNIFORM_TEMPERATURE: f64 = 1e-3;




/**
 * Explosion in a horizontally stratified external medium
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HaloKilonova {
    pub altitude: f64,
    pub external_medium_density: f64,
    pub launch_radius: f64,
    pub shell_thickness: f64,
    pub kinetic_energy: f64,
    pub shell_mass: f64,
    pub radial_distance: f64,
}


// ============================================================================
impl HaloKilonova {


    /**
     * Return the radial extent (in cm) of the shell at time t.
     */
    fn shell_extent(&self, t: f64) -> std::ops::Range<f64> {
        let r_outer_shell_surface = self.launch_radius + self.shell_velocity() * t;
        let r_inner_shell_surface = self.launch_radius + self.shell_velocity() * (t - self.shell_duration());
        r_inner_shell_surface..r_outer_shell_surface
    }


    /**
     * The velocity (in cm/s) the shell moves at (computed from the shell mass
     * and kinetic energy). Note this is expression assumes the shell is
     * sub-relativistic.
     */
    fn shell_velocity(&self) -> f64 {
        (2.0 * self.kinetic_energy / self.shell_mass).sqrt()
    }


    /**
     * The duration (in s) during which the shell is emerging from the inner
     * boundary.
     */
    fn shell_duration(&self) -> f64 {
        self.shell_thickness / self.shell_velocity()
    }

}

// ============================================================================
impl InitialModel for HaloKilonova {

    fn validate(&self) -> anyhow::Result<()> {
        if self.shell_velocity() > 0.25 * LIGHT_SPEED {
            anyhow::bail!{"
             The shell is moving faster (v/c = {}) than 0.25 c, but
             this problem assumes Newtonian expressions for the
             kinetic energy. Consider reducing the kinetic energy or
             increasing the shell mass.", self.shell_velocity() / LIGHT_SPEED}
        } else {
            Ok(())
        }
    }

    fn primitive_at(&self, coordinate: (f64, f64), t: f64) -> AnyPrimitive {
        
        let (r, q) = coordinate;

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
        } else {
            let z = r * q.cos() + self.altitude;
            let d0 = self.external_medium_density;
            let d   =  if z < 0.0 {
                d0
            } else {
                GalacticModel::density(&GalacticModel{g: 6.67e-8, m_b: 3.377e43,
                                                      a_b: 8.98e20, v_h: 1.923e7, 
                                                      a_h: 9.26e22, m_s: 1.538e44,
                                                      a_s: 1.461e22, b_s: 1.790e21,
                                                      m_g: 5.434e43, a_g: 1.461e22,
                                                      b_g: 7.035e23}, self.radial_distance, z).thin_disk
            };
            let p   = d*1.0e-3;

            AnyPrimitive {
                velocity_r: 0.0,
                velocity_q: 0.0,
                mass_density: d,
                gas_pressure: p,
            }            
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
