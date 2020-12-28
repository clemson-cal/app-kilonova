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
pub trait Conserved: 'static + Clone + Copy + Send + Sync + Arithmetic {
    fn lab_frame_mass(&self) -> f64;
}


/**
 * Primitive field type for the hydrodynamics system
 */
pub trait Primitive: 'static + Clone + Copy + Send + Sync + Arithmetic + Serialize + Default {
} 


/**
 * Interface to a hydrodynamics system: either euler_2d or srhd_2d
 */
pub trait Hydrodynamics: Clone {

    /// The type of the conserved struct: mass, momentum, energy
    type Conserved: Conserved;

    /// The type of the primitive struct: comoving density, (four)-velocity, pressure
    type Primitive: Primitive;

    /// Compute the PLM difference from a stencil of colinear primitive states
    fn plm_gradient(&self, theta: f64, a: &Self::Primitive, b: &Self::Primitive, c: &Self::Primitive) -> Self::Primitive;

    /// Convert from a primitive to a conserved state
    fn to_primitive(&self, u: Self::Conserved) -> Self::Primitive;

    /// Convert from a conserved to a primitive state (signature may be
    /// changed to return Result)
    fn to_conserved(&self, p: Self::Primitive) -> Self::Conserved;

    /// Convert from an agnostic primitive state to the one specific to this
    /// hydrodynamics system
    fn interpret(&self, agnostic: &AgnosticPrimitive) -> Self::Primitive;


    fn intercell_flux(&self, pl: Self::Primitive, pr: Self::Primitive, sl: f64, sr: f64, direction: Direction) -> Self::Conserved;
}


/**
 * Implemented by types that can generate primitive fields to be used as an
 * initial or boundary value
 */
#[enum_dispatch::enum_dispatch]
pub trait InitialModel: Clone {

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
