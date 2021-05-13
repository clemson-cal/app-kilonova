use serde::{Serialize, Deserialize};
use godunov_core::piecewise_linear;
use godunov_core::runge_kutta::RungeKuttaOrder;
use crate::physics::{AnyPrimitive, RiemannSolver, Direction, HydroErrorType, LIGHT_SPEED};
use crate::traits::Hydrodynamics;
use crate::galmod::GalacticModel;




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

    /// Define the time step based on the maximum signal speed. If false,
    /// assume the speed of light.
    #[serde(default)]
    pub adaptive_time_step: bool,
}




// ============================================================================
impl Hydrodynamics for RelativisticHydro {
    type Conserved = hydro_srhd::srhd_2d::Conserved;
    type Primitive = hydro_srhd::srhd_2d::Primitive;

    fn validate(&self) -> anyhow::Result<()> {
        if self.plm_theta < 1.0 || self.plm_theta > 2.0 {
            anyhow::bail!("plm_theta must be in the range [1, 2]")            
        }
        // commented out for now to allow to force a particular time step
        // if self.cfl_number < 0.0 || self.cfl_number > 0.7 {
        //     anyhow::bail!("cfl_number must be in the range [0.0, 0.7]")
        // }
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

    fn try_to_primitive(&self, u:Self::Conserved) -> Result<Self::Primitive, HydroErrorType>{

        if u.lab_frame_density() < 0.0 {
            return Err(HydroErrorType::NegativeDensity(u.lab_frame_density()))
        }
        else if u.energy_density() < 0.0 {
            return Err(HydroErrorType::NegativeEnergyDensity(u.energy_density()))
        }

        let valid_primitive = match u.to_primitive(self.gamma_law_index) {
            hydro_srhd::srhd_2d::RecoveredPrimitive::Success(p) => p,
            hydro_srhd::srhd_2d::RecoveredPrimitive::NegativePressure(p) => {
                hydro_srhd::srhd_2d::Primitive(p.0, p.1, p.2, 1e-3 * p.0)
            }
            hydro_srhd::srhd_2d::RecoveredPrimitive::RootFinderFailed(u) => {
                return Err(HydroErrorType::RootFinderFailed(u))?
            }
        };

        Ok(valid_primitive)
    }

    fn to_primitive(&self, u: Self::Conserved) -> Self::Primitive {
        self.try_to_primitive(u).unwrap()
    }

    fn to_conserved(&self, p: Self::Primitive) -> Self::Conserved {
        p.to_conserved(self.gamma_law_index)
    }

    fn max_signal_speed(&self, p: Self::Primitive) -> f64 {
        p.max_signal_speed(self.gamma_law_index) * LIGHT_SPEED
    }

    fn global_signal_speed(&self) -> Option<f64> {
        if self.adaptive_time_step {
            None
        } else {
            Some(LIGHT_SPEED)
        }
    }

    fn interpret(&self, a: &AnyPrimitive) -> Self::Primitive {
        hydro_srhd::srhd_2d::Primitive(a.mass_density, a.velocity_r, a.velocity_q, a.gas_pressure)
    }

    fn any(&self, p: &Self::Primitive) -> AnyPrimitive {
        AnyPrimitive {
            velocity_r: p.gamma_beta_1(),
            velocity_q: p.gamma_beta_2(),
            mass_density: p.mass_density(),
            gas_pressure: p.gas_pressure(),
        }
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

    fn gravitational_source_terms(&self, p: Self::Primitive, coordinate: (f64, f64)) -> Self::Conserved {
        let h0 = p.specific_enthalpy(self.gamma_law_index);
        let gmod = GalacticModel {g: 6.67e-8,
                                   m_b: 3.377e43,
                                   a_b: 8.98e20,
                                   v_h: 1.923e7,
                                   a_h: 9.26e22,
                                   m_s: 1.538e44,
                                   a_s: 1.461e22,
                                   b_s: 1.790e21,
                                   m_g: 5.434e43,
                                   a_g: 1.461e22,
                                   b_g: 7.035e23,
                                  };
        let cosq = f64::cos(coordinate.1);
        let sinq = f64::sin(coordinate.1);
        let gz = gmod.g_field_z(1e22, coordinate.0*cosq + 1.5e20).total();
        
        let gd = 0.0;
        let gr = p.lorentz_factor() * p.mass_density() * h0 * gz * cosq / LIGHT_SPEED / LIGHT_SPEED;
        let gq = -p.lorentz_factor() * p.mass_density() * h0 * gz * sinq / LIGHT_SPEED / LIGHT_SPEED;
        let ge = p.lorentz_factor() * p.mass_density() * h0 * gz * cosq / LIGHT_SPEED * (p.gamma_beta_1()*cosq - p.gamma_beta_2()*sinq);

        hydro_srhd::srhd_2d::Conserved(gd, gr, gq, ge)
    }

    fn cfl_number(&self) -> f64 {
        self.cfl_number
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
    fn lorentz_factor(&self) -> f64 {
        self.lorentz_factor_squared().sqrt()
    }
}
