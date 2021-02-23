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
pub struct AnyPrimitive {

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
impl Into<[f64; 4]> for AnyPrimitive {
    fn into(self) -> [f64; 4] {
        [self.velocity_r, self.velocity_q, self.mass_density, self.gas_pressure]
    }
}

impl From<[f64; 4]> for AnyPrimitive {
    fn from(d: [f64; 4]) -> Self {
        AnyPrimitive{
            velocity_r: d[0],
            velocity_q: d[1],
            mass_density: d[2],
            gas_pressure: d[3],
        }
    }
}

/// Error Implementation
// ============================================================================
type Conserved = hydro_srhd::srhd_2d::Conserved;

#[derive(thiserror::Error, Debug, Clone)]
pub enum HydroErrorType {
    #[error("Negative Mass Density {0:.4e}")]
    NegativeDensity(f64),

    #[error("Negative Pressure {0:.4e}")]
    NegativePressure(f64),

    #[error("The Root Finder Failed to Converge {0:?}")]
    RootFinderFailed(Conserved),
}

impl HydroErrorType {
    pub fn at_position(self, position: (f64, f64)) -> HydroError {
        HydroError{source: self, position}
    }
}


#[derive(thiserror::Error, Debug, Clone)]
#[error("at position ({:.4} {:.4}) in the mesh",
    position.0,
    position.1,
)]

// ============================================================================
pub struct HydroError {
    source: HydroErrorType,
    position: (f64, f64),
}

impl HydroError {
    pub fn with_model(self) -> Self {
        Self {
            source: self.source,
            position: self.position,
        }
    }
}