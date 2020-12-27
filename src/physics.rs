use serde::{Serialize, Deserialize};
use godunov_core::piecewise_linear;
use crate::traits::Hydrodynamics;




/**
 * Primitive variable state that is agnostic to the hydrodynamics system
 */
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




/**
 * Interface implementation for relativistic hydrodynamics
 */
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, tag = "type")]
pub struct RelativisticHydrodynamics {
    pub gamma_law_index: f64,
}




// ============================================================================
impl Hydrodynamics for RelativisticHydrodynamics {
    type Conserved = hydro_srhd::srhd_2d::Conserved;
    type Primitive = hydro_srhd::srhd_2d::Primitive;

    fn plm_gradient(&self, theta: f64, a: &Self::Primitive, b: &Self::Primitive, c: &Self::Primitive) -> Self::Primitive {
        piecewise_linear::plm_gradient4(theta, a, b, c)
    }

    fn to_primitive(&self, u: Self::Conserved) -> Self::Primitive {
        u.to_primitive(self.gamma_law_index).unwrap()
    }

    fn to_conserved(&self, p: Self::Primitive) -> Self::Conserved {
        p.to_conserved(self.gamma_law_index)
    }

    fn interpret(&self, a: &AgnosticPrimitive) -> Self::Primitive {
        hydro_srhd::srhd_2d::Primitive(a.mass_density, a.velocity_r, a.velocity_q, a.gas_pressure)
    }
}




// ============================================================================
impl crate::traits::Arithmetic for hydro_srhd::srhd_2d::Conserved {
}

impl crate::traits::Conserved for hydro_srhd::srhd_2d::Conserved {
    fn lab_frame_mass(&self) -> f64 {
        self.lab_frame_density()
    }
}

impl crate::traits::Primitive for hydro_srhd::srhd_2d::Primitive {
}
