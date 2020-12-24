use std::ops::{Add, Sub, Mul, Div};




// ============================================================================
pub trait Arithmetic: Add<Output=Self> + Sub<Output=Self> + Mul<f64, Output=Self> + Div<f64, Output=Self> + Sized {}




// ============================================================================
pub trait Zeros {
    fn zeros() -> Self;
}




// ============================================================================
pub trait Conserved: 'static + Clone + Copy + Send + Sync + Zeros + Arithmetic { //+ hdf5::H5Type {
    fn mass_and_momentum(&self) -> (f64, f64, f64);
}




// ============================================================================
pub trait Primitive: Clone + Copy + Send + Sync { //+ hdf5::H5Type {
    fn velocity_x(self) -> f64;
    fn velocity_y(self) -> f64;
    fn mass_density(self) -> f64;
}




// ============================================================================
pub trait Hydrodynamics: Copy + Send
{
    type Conserved: Conserved;
    type Primitive: Primitive;

    fn gamma_law_index(&self) -> f64;
    fn plm_gradient(&self, theta: f64, a: &Self::Primitive, b: &Self::Primitive, c: &Self::Primitive) -> Self::Primitive;
    fn to_primitive(&self, u: Self::Conserved) -> Self::Primitive;
    fn to_conserved(&self, p: Self::Primitive) -> Self::Conserved;
}
