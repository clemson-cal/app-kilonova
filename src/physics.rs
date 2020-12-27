use ndarray::{Array, Ix2};
use godunov_core::piecewise_linear;
use crate::traits::{Hydrodynamics, InitialModel, Primitive};
use crate::mesh::SphericalPolarGrid;




/**
 * Primitive state that is agnostic to the hydrodynamics system
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
pub struct RelativisticHydrodynamics {
    pub gamma_law_index: f64,
}




// ============================================================================
// impl RelativisticHydrodynamics {
//     pub fn new() -> Self {
//         Self{
//             gamma_law_index: 4.0 / 3.0
//         }
//     }
// }


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
}

impl crate::traits::Primitive for hydro_srhd::srhd_2d::Primitive {
}




/**
 * Generate an array of hydrodynamic primitives for the given grid patch and
 * model.
 */
pub fn grid_primitive<H, M, P>(grid: &SphericalPolarGrid, system: &H, model: &M) -> Array<P, Ix2>
where
    H: Hydrodynamics<Primitive = P>,
    P: Primitive,
    M: InitialModel {
    Array::from_shape_fn(grid.cell_dim(), |(i, j)| {
        system.interpret(&model.at(grid.zone(i as i64, j as i64).centroid()))
    })
}
