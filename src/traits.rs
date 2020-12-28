use std::ops::{Add, Sub, Mul, Div};
use serde::Serialize;
use crate::physics::AgnosticPrimitive;


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
pub trait Primitive: 'static + Clone + Copy + Send + Sync + Serialize {
} 


/**
 * Interface to a hydrodynamics system: either euler_2d or srhd_2d
 */
pub trait Hydrodynamics: Clone {
    type Conserved: Conserved;
    type Primitive: Primitive;

    fn plm_gradient(&self, theta: f64, a: &Self::Primitive, b: &Self::Primitive, c: &Self::Primitive) -> Self::Primitive;
    fn to_primitive(&self, u: Self::Conserved) -> Self::Primitive;
    fn to_conserved(&self, p: Self::Primitive) -> Self::Conserved;
    fn interpret(&self, agnostic: &AgnosticPrimitive) -> Self::Primitive;
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
     fn primitive_at(&self, coordinate: (f64, f64)) -> AgnosticPrimitive;

     /**
      * Return the scalar concentration at the given r-theta coordinate.
      */
     fn scalar_at(&self, coordinate: (f64, f64)) -> f64;
}
