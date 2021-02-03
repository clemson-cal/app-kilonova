mod relativistic_hydro;
mod newtonian_hydro;

use serde::{Serialize, Deserialize};
pub use relativistic_hydro::RelativisticHydro;
pub use newtonian_hydro::NewtonianHydro;
pub static LIGHT_SPEED: f64 = 3e10;




/**
 * Enum for the cardinal grid axes
 */
pub enum Direction {
    Polar,
    Radial,
}




/**
 * Enum for Riemann solver type
 */
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub enum RiemannSolver {
    HLLE,
    HLLC,
}




/**
 * Primitive variable state that is agnostic to the hydrodynamics system
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(from = "[f64; 4]", into = "[f64; 4]")]
pub struct AgnosticPrimitive {

    /// Radial velocity (radial gamma-beta for relativistic hydro)
    pub velocity_r: f64,

    /// Polar velocity (polar gamma-beta for relativistic hydro)
    pub velocity_q: f64,

    /// Mass density (comoving for relativistic)
    pub mass_density: f64,

    /// Gas pressure
    pub gas_pressure: f64,
}




// ============================================================================
impl Into<[f64; 4]> for AgnosticPrimitive {
    fn into(self) -> [f64; 4] {
        [self.velocity_r, self.velocity_q, self.mass_density, self.gas_pressure]
    }
}

impl From<[f64; 4]> for AgnosticPrimitive {
    fn from(d: [f64; 4]) -> Self {
        AgnosticPrimitive{
            velocity_r: d[0],
            velocity_q: d[1],
            mass_density: d[2],
            gas_pressure: d[3],
        }
    }
}
