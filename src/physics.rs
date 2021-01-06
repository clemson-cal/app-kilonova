use serde::{Serialize, Deserialize};
use godunov_core::piecewise_linear;
use godunov_core::runge_kutta::RungeKuttaOrder;
use crate::mesh::Mesh;
use crate::state::State;
use crate::traits::Hydrodynamics;

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
#[derive(Clone, Serialize, Deserialize)]
pub enum RiemannSolver {
    HLLE,
    HLLC,
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

    /// Index for the gamma-law equation of state
    pub gamma_law_index: f64,

    /// Parameter for gradient estimation: [1, 2]
    pub plm_theta: f64,

    /// Time step size: [0.0, 0.7]
    pub cfl_number: f64,

    /// Runge-Kutta order: [RK1 | RK2 | RK3]
    pub runge_kutta_order: RungeKuttaOrder,

    /// Riemann solver: [HLLE | HLLC]
    pub riemann_solver: RiemannSolver,
}




// ============================================================================
impl Hydrodynamics for RelativisticHydro {
    type Conserved = hydro_srhd::srhd_2d::Conserved;
    type Primitive = hydro_srhd::srhd_2d::Primitive;

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

    fn time_step(&self, state: &State<Self::Conserved>, mesh: &Mesh) -> f64 {
        let (index, ..) = state.inner_outer_block_indexes();
        self.cfl_number * mesh.smallest_spacing(index) / LIGHT_SPEED
    }

    fn plm_gradient_primitive(&self, a: &Self::Primitive, b: &Self::Primitive, c: &Self::Primitive) -> Self::Primitive {
        piecewise_linear::plm_gradient4(self.plm_theta, a, b, c)
    }

    fn plm_gradient_scalar(&self, a: &f64, b: &f64, c: &f64) -> f64 {
        piecewise_linear::plm_gradient(self.plm_theta, a, b, c)
    }

    fn to_primitive(&self, u: Self::Conserved) -> Self::Primitive {
        match u.to_primitive(self.gamma_law_index) {
            hydro_srhd::srhd_2d::RecoveredPrimitive::Success(p) => p,
            hydro_srhd::srhd_2d::RecoveredPrimitive::NegativePressure(p) => {
                hydro_srhd::srhd_2d::Primitive(p.0, p.1, p.2, 1e-3 * p.0)
            },
            hydro_srhd::srhd_2d::RecoveredPrimitive::RootFinderFailed(u) => panic!("c2p root finder failed {:?}", u),
        }
    }

    fn to_conserved(&self, p: Self::Primitive) -> Self::Conserved {
        p.to_conserved(self.gamma_law_index)
    }

    fn interpret(&self, a: &AgnosticPrimitive) -> Self::Primitive {
        hydro_srhd::srhd_2d::Primitive(a.mass_density, a.velocity_r, a.velocity_q, a.gas_pressure)
    }

    fn intercell_flux(&self, pl: Self::Primitive, pr: Self::Primitive, sl: f64, sr: f64, direction: Direction) -> (Self::Conserved, f64) {
        let mode = match self.riemann_solver {
            RiemannSolver::HLLE => hydro_srhd::srhd_2d::RiemannSolverMode::HlleFlux,
            RiemannSolver::HLLC => hydro_srhd::srhd_2d::RiemannSolverMode::HllcFlux,
        };            
        let axis = match direction {
            Direction::Radial => hydro_srhd::geometry::Direction::X,
            Direction::Polar  => hydro_srhd::geometry::Direction::Y,
        };
        let (f, g, _) = hydro_srhd::srhd_2d::riemann_hllc_scalar(pl, pr, sl, sr, axis, self.gamma_law_index, mode);
        (f * LIGHT_SPEED, g * LIGHT_SPEED)
    }

    fn geometrical_source_terms(&self, p: Self::Primitive, coordinate: (f64, f64)) -> Self::Conserved {
        p.spherical_geometry_source_terms(coordinate.0, coordinate.1, self.gamma_law_index) * LIGHT_SPEED
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
