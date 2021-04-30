use serde::{Serialize, Deserialize};
use godunov_core::piecewise_linear;
use godunov_core::runge_kutta::RungeKuttaOrder;
use crate::physics::{AnyPrimitive, Direction, HydroErrorType};
use crate::traits::Hydrodynamics;




/**
 * Interface implementation for Newtonian hydrodynamics
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NewtonianHydro {

    /// Index for the gamma-law equation of state
    pub gamma_law_index: f64,

    /// Parameter for gradient estimation: [1, 2]
    pub plm_theta: f64,

    /// Time step size: [0.0, 0.7]
    pub cfl_number: f64,

    /// Runge-Kutta order: [RK1 | RK2 | RK3]
    pub runge_kutta_order: RungeKuttaOrder,
}




// ============================================================================
impl Hydrodynamics for NewtonianHydro {
    type Conserved = hydro_euler::euler_2d::Conserved;
    type Primitive = hydro_euler::euler_2d::Primitive;

    fn validate(&self) -> anyhow::Result<()> {
        if self.plm_theta < 1.0 || self.plm_theta > 2.0 {
            anyhow::bail!("plm_theta must be in the range [1, 2]")            
        }
        if self.cfl_number < 0.0 || self.cfl_number > 0.7 {
            anyhow::bail!("cfl_number must be in the range [0.0, 0.7]")
        }
        Ok(())
    }

    fn runge_kutta_order(&self) -> RungeKuttaOrder {
        self.runge_kutta_order
    }

    fn plm_gradient_primitive(&self, a: &Self::Primitive, b: &Self::Primitive, c: &Self::Primitive) -> Self::Primitive {
        piecewise_linear::plm_gradient4(self.plm_theta, a, b, c)
    }

    fn plm_gradient_scalar(&self, a: &f64, b: &f64, c: &f64) -> f64 {
        piecewise_linear::plm_gradient(self.plm_theta, a, b, c)
    }

    fn try_to_primitive(&self, u: Self::Conserved) -> Result<Self::Primitive, HydroErrorType> {
        if u.mass_density() < 0.0 {
            return Err(HydroErrorType::NegativeDensity(u.mass_density()))
        }
        Ok(u.to_primitive(self.gamma_law_index))
    }

    fn to_primitive(&self, u: Self::Conserved) -> Self::Primitive {
        self.try_to_primitive(u).unwrap()
    }

    fn to_conserved(&self, p: Self::Primitive) -> Self::Conserved {
        p.to_conserved(self.gamma_law_index)
    }

    fn max_signal_speed(&self, p: Self::Primitive) -> f64 {
        p.max_signal_speed(self.gamma_law_index)
    }

    fn global_signal_speed(&self) -> Option<f64> {
        None
    }

    fn interpret(&self, a: &AnyPrimitive) -> Self::Primitive {
        hydro_euler::euler_2d::Primitive(a.mass_density, a.velocity_r, a.velocity_q, a.gas_pressure)
    }

    fn any(&self, p: &Self::Primitive) -> AnyPrimitive {
        AnyPrimitive {
            velocity_r: p.velocity_1(),
            velocity_q: p.velocity_2(),
            mass_density: p.mass_density(),
            gas_pressure: p.gas_pressure(),
        }
    }

    fn intercell_flux(&self, pl: Self::Primitive, pr: Self::Primitive, sl: f64, sr: f64, direction: Direction) -> (Self::Conserved, f64) {
        let axis = match direction {
            Direction::Radial => hydro_euler::geometry::Direction::X,
            Direction::Polar  => hydro_euler::geometry::Direction::Y,
        };
        hydro_euler::euler_2d::riemann_hlle_scalar(pl, pr, sl, sr, axis, self.gamma_law_index)
    }

    fn geometrical_source_terms(&self, p: Self::Primitive, coordinate: (f64, f64)) -> Self::Conserved {
        p.spherical_geometry_source_terms(coordinate.0, coordinate.1)
    }

    fn gravitational_source_terms(&self, _p: Self::Primitive, _coordinate: (f64, f64)) -> Self::Conserved {
        let gd = 0.0;
        let gr = 0.0;
        let gq = 0.0;
        let ge = 0.0;
        hydro_euler::euler_2d::Conserved(gd, gr, gq, ge)
    }

    fn cfl_number(&self) -> f64 {
        self.cfl_number
    }
}




// ============================================================================
impl crate::traits::Arithmetic for hydro_euler::euler_2d::Conserved {
}

impl crate::traits::Conserved for hydro_euler::euler_2d::Conserved {
    fn lab_frame_mass(&self) -> f64 {
        self.mass_density()
    }
}

impl crate::traits::Arithmetic for hydro_euler::euler_2d::Primitive {
}

impl crate::traits::Primitive for hydro_euler::euler_2d::Primitive {
    fn lorentz_factor(&self) -> f64 {
        1.0
    }
}
