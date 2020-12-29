use std::ops::{Add, Sub, Mul, Div};
use serde::Serialize;
use crate::physics::{AgnosticPrimitive, Direction};


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
} 


/**
 * Interface to a hydrodynamics system: either euler_2d or srhd_2d
 */
pub trait Hydrodynamics: Clone {

    /// The type of the conserved struct: mass, momentum, energy.
    type Conserved: Conserved;

    /// The type of the primitive struct: comoving density, (four)-velocity,
    /// pressure.
    type Primitive: Primitive;

    /// Compute the PLM difference from a stencil of colinear primitive
    /// states
    fn plm_gradient_primitive(&self, a: &Self::Primitive, b: &Self::Primitive, c: &Self::Primitive) -> Self::Primitive;

    /// Compute the PLM difference from a stencil of colinear scalar
    /// concentration states
    fn plm_gradient_scalar(&self, a: &f64, b: &f64, c: &f64) -> f64;

    /// Convert from a primitive to a conserved state.
    fn to_primitive(&self, u: Self::Conserved) -> Self::Primitive;

    /// Convert from a conserved to a primitive state (signature may be
    /// changed to return Result).
    fn to_conserved(&self, p: Self::Primitive) -> Self::Conserved;

    /// Convert from an agnostic primitive state to the one specific to this
    /// hydrodynamics system.
    fn interpret(&self, agnostic: &AgnosticPrimitive) -> Self::Primitive;

    /// Return the Godunov flux of the conserved quantities and the passive
    /// scalar, given reconstructued values of the primitives (`pl`, `pr`) and
    /// the scalar concentration (`sl`, `sr`) to either side of the cell
    /// interface.
    fn intercell_flux(&self, pl: Self::Primitive, pr: Self::Primitive, sl: f64, sr: f64, direction: Direction) -> (Self::Conserved, f64);

    /// Return the geometrical source terms (conserved quantity per unit
    /// volume) for the given primitive state and r-theta coordinate.
    fn geometrical_source_terms(&self, p: Self::Primitive, coordinate: (f64, f64)) -> Self::Conserved;
}


/**
 * Implemented by types that can generate primitive fields to be used as an
 * initial or boundary value
 */
#[enum_dispatch::enum_dispatch]
pub trait InitialModel: Clone {

    /**
     * Return an error if this model was not configured to yield acceptable
     * data.
     */
    fn validate(&self) -> anyhow::Result<()>;

    /**
     * Return an agnostic primitive state at the give r-theta coordinate. An
     * [`AgnosticPrimitive`] must be converted to the appropriate [`Primitive`]
     * type by the [`Hydrodynamics::interpret`] method.
     */
     fn primitive_at(&self, coordinate: (f64, f64), time: f64) -> AgnosticPrimitive;

     /**
      * Return the scalar concentration at the given r-theta coordinate.
      */
     fn scalar_at(&self, coordinate: (f64, f64), time: f64) -> f64;
}
