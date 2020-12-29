use serde::{Serialize, Deserialize};
use godunov_core::piecewise_linear;
use crate::traits::Hydrodynamics;




pub enum Direction {
    Polar,
    Radial,
}




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
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RelativisticHydro {
    pub gamma_law_index: f64,
    pub plm_theta: f64,
}




// ============================================================================
impl Hydrodynamics for RelativisticHydro {
    type Conserved = hydro_srhd::srhd_2d::Conserved;
    type Primitive = hydro_srhd::srhd_2d::Primitive;

    fn plm_gradient_primitive(&self, a: &Self::Primitive, b: &Self::Primitive, c: &Self::Primitive) -> Self::Primitive {
        piecewise_linear::plm_gradient4(self.plm_theta, a, b, c)
    }

    fn plm_gradient_scalar(&self, a: &f64, b: &f64, c: &f64) -> f64 {
        piecewise_linear::plm_gradient(self.plm_theta, a, b, c)
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

    fn intercell_flux(&self, pl: Self::Primitive, pr: Self::Primitive, sl: f64, sr: f64, direction: Direction) -> (Self::Conserved, f64) {
        let mode = hydro_srhd::srhd_2d::RiemannSolverMode::HlleFlux;
        let axis = match direction {
            Direction::Radial => hydro_srhd::geometry::Direction::X,
            Direction::Polar  => hydro_srhd::geometry::Direction::Y,
        };
        let (f, g, _) = hydro_srhd::srhd_2d::riemann_hllc_scalar(pl, pr, sl, sr, axis, self.gamma_law_index, mode);
        (f, g)
    }

    fn geometrical_source_terms(&self, p: Self::Primitive, coordinate: (f64, f64)) -> Self::Conserved {
        p.spherical_geometry_source_terms(coordinate.0, coordinate.1, self.gamma_law_index)
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

impl crate::traits::Arithmetic for hydro_srhd::srhd_2d::Primitive {
}

impl crate::traits::Primitive for hydro_srhd::srhd_2d::Primitive {
}
