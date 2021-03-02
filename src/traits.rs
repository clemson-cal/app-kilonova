use std::ops::{Add, Sub, Mul, Div};
use serde::Serialize;
use godunov_core::runge_kutta::RungeKuttaOrder;
use crate::mesh::Mesh;
use crate::state::State;
use crate::physics::{AnyPrimitive, Direction, HydroError, HydroErrorType};




/**
 * Implemented by types vector-like types
 */
pub trait Arithmetic: Add<Output=Self> + Sub<Output=Self> + Mul<f64, Output=Self> + Div<f64, Output=Self> + Sized {
}




/**
 * Conserved field type for the hydrodynamics system
 */
pub trait Conserved: 'static + Clone + Copy + Send + Sync + Arithmetic + Default {
    fn lab_frame_mass(&self) -> f64;
}




/**
 * Primitive field type for the hydrodynamics system
 */
pub trait Primitive: 'static + Clone + Copy + Send + Sync + Arithmetic + Default + Serialize {
    fn lorentz_factor(&self) -> f64;
}




/**
 * Interface to a hydrodynamics system: either euler_2d or srhd_2d
 */
pub trait Hydrodynamics: 'static + Clone + Send {

    /**
     * The type of the conserved struct: mass, momentum, energy.
     */
    type Conserved: Conserved;

    /**
     * The type of the primitive struct: comoving density, (four)-velocity,
     * pressure.
     */
    type Primitive: Primitive;

    /**
     * Return an error if the hydrodynamics instance was somehow configured
     * improperly.
     */
    fn validate(&self) -> anyhow::Result<()>;

    /**
     * Return the Runge Kutta order, which may be a user configurable value.
     */
    fn runge_kutta_order(&self) -> RungeKuttaOrder;

    /**
     * Compute the PLM difference from a stencil of colinear primitive
     * states.
     */
    fn plm_gradient_primitive(&self, a: &Self::Primitive, b: &Self::Primitive, c: &Self::Primitive) -> Self::Primitive;

    /**
     * Compute the PLM difference from a stencil of colinear scalar
     * concentration states.
     */
    fn plm_gradient_scalar(&self, a: &f64, b: &f64, c: &f64) -> f64;

    /**
     * Try to convert from a conserved to a primitive hydrodynamic state,
     * returning an appropriate error type if the conversion failed. This
     * function is not permitted to panic.
     */
    fn try_to_primitive(&self, u: Self::Conserved) -> Result<Self::Primitive, HydroErrorType>;

    /**
     * Convert from a conserved to a primitive hydrodynamic state. This function
     * is is permitted to panic if the conversion fails.
     */
    fn to_primitive(&self, u: Self::Conserved) -> Self::Primitive;

    /**
     * Convert from a primitive to a conserved state.
     */
    fn to_conserved(&self, p: Self::Primitive) -> Self::Conserved;

    fn max_signal_speed(&self, p: Self::Primitive) -> f64;

    /**
     * Convert from an any-primitive state to the one specific to this
     * hydrodynamics system.
     */
    fn interpret(&self, any: &AnyPrimitive) -> Self::Primitive;

    /**
     * Convert from a primitive state specific to this hydrodynamics system to
     * an any-primitive.
     */
    fn any(&self, p: &Self::Primitive) -> AnyPrimitive;

    /**
     * Return the Godunov flux of the conserved quantities and the passive
     * scalar, given reconstructued values of the primitives (`pl`, `pr`) and
     * the scalar concentration (`sl`, `sr`) to either side of the cell
     * interface.
     */
    fn intercell_flux(&self, pl: Self::Primitive, pr: Self::Primitive, sl: f64, sr: f64, direction: Direction) -> (Self::Conserved, f64);

    /**
     * Return the geometrical source terms (conserved quantity per unit volume)
     * for the given primitive state and r-theta coordinate.
     */
    fn geometrical_source_terms(&self, p: Self::Primitive, coordinate: (f64, f64)) -> Self::Conserved;

    /**
     * Return the CFL number to be used
     */
    fn cfl_number(&self) -> f64;

    /**
     * Return the time step size, computed from the mesh, the hydrodynamics
     * state, and internal parameters such as the CFL number.
     */
    fn time_step(&self, state: &State<Self::Conserved>, mesh: &Mesh) -> Result<f64, HydroError> {
        Ok(state.solution.iter().try_fold(f64::MAX, |dt, (index, state)| {
            let geometry = mesh.subgrid(*index).geometry();
            let block_dt = state
                .try_to_primitive(self, &geometry)?
                .iter()
                .zip(&geometry.cell_linear_dimension())
                .fold(dt, |dt, (p, dl)| dt.min(dl / self.max_signal_speed(*p))
            );
            Ok(dt.min(block_dt))
        })? * self.cfl_number())
    }
}




/**
 * Implemented by types that can generate primitive fields to be used as an
 * initial or boundary value
 */
pub trait InitialModel: Clone {

    /**
     * Return an error if this model was not configured to yield acceptable
     * data.
     */
    fn validate(&self) -> anyhow::Result<()>;

    /**
     * Return an agnostic primitive state at the give r-theta coordinate. An
     * [`AnyPrimitive`] must be converted to the appropriate [`Primitive`]
     * type by the [`Hydrodynamics::interpret`] method.
     */
     fn primitive_at(&self, coordinate: (f64, f64), time: f64) -> AnyPrimitive;

     /**
      * Return the scalar concentration at the given r-theta coordinate.
      */
     fn scalar_at(&self, coordinate: (f64, f64), time: f64) -> f64;
}
