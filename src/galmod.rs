//! 10-parameter axisymmetric model of a disk galaxy, using
//! the Plummer model for the central bulge, two Miyamoto-
//! Nagai models for the thin and thick disks, and a logarithmic
//! model for the dark matter halo.

use std::f64::consts::PI;
use derive_more::{Add,Sub, Mul, Div};

// Define struct for galactic model parameters,
// including the gravitational constant. Here, slr
// stands for solar masses.
#[derive(Debug, Add, Sub, Mul, Div)]
pub struct GalacticModel {
    pub g:   f64, // gravitational constant (kpc*kpc*kpc/Myr/Myr/slr)
    pub m_b: f64, // mass of central bulge (slr)
    pub a_b: f64, // radial scale length of central bulge (kpc)
    pub v_h: f64, // radial velocities at large distances (kpc/Myr)
    pub a_h: f64, // radial scale length of dark matter halo (kpc)
    pub m_s: f64, // mass of thin disk (slr)
    pub a_s: f64, // radial scale length of thin disk (kpc)
    pub b_s: f64, // vertical scale length of thin disk (kpc)
    pub m_g: f64, // mass of thick disk (slr)
    pub a_g: f64, // radial scale length of thick disk (kpc)
    pub b_g: f64, // vertical scale length of thick disk (kpc)
}

impl GalacticModel {

    // Gravitational potential
    pub fn potential(&self, r: f64, z: f64) -> f64 {

        let phi_bulge      = -self.g * self.m_b / (r*r + z*z + self.a_b*self.a_b).sqrt();
        let phi_thin_disk  = -self.g * self.m_s / (r*r + (self.a_s + (z*z + self.b_s*self.b_s).sqrt()).powi(2)).sqrt();
        let phi_thick_disk = -self.g * self.m_g / (r*r + (self.a_g + (z*z + self.b_g*self.b_g).sqrt()).powi(2)).sqrt();
        let phi_halo       = 0.5*(self.v_h).powi(2)*((r*r + z*z + self.a_h*self.a_h).ln());

        phi_bulge + phi_thin_disk + phi_thick_disk + phi_halo

    }

    // Mass density profile obtained via Poisson's equation using the above potential
    pub fn density(&self, r: f64, z: f64) -> f64 {

        let rho_bulge        = (3.0*self.m_b/4.0/PI/(self.a_b).powi(3))*
                               (1.0 + (r*r + z*z)/self.a_b/self.a_b).powf(-5.0/2.0);

        let rho_thin_disk_z  = -self.m_s/4.0/PI*((4.0*self.a_s*z.powi(4) + (z*z + self.b_s*self.b_s).sqrt()
                               *((2.0*self.a_s*self.a_s - r*r)*z*z - self.b_s*self.b_s*r*r - 3.0*self.a_s*self.a_s*self.b_s*self.b_s)
                               + (z*z + self.b_s*self.b_s).powf(3.0/2.0)*(2.0*z*z - self.b_s*self.b_s) + self.a_s*self.b_s*self.b_s*z*z
                               - self.a_s*self.b_s*self.b_s*r*r - 3.0*self.a_s*(self.b_s).powi(4)-(self.a_s).powi(3)*self.b_s*self.b_s)
                               /((z*z + self.b_s*self.b_s).powf(3.0/2.0)*(((z*z + self.b_s*self.b_s).sqrt() + self.a_s).powi(2) + r*r).powf(5.0/2.0)));
        let rho_thin_disk_r  = self.m_s/4.0/PI*((2.0*(r*r + (self.a_s + (z*z + self.b_s*self.b_s).sqrt()).powi(2)) - 3.0*r*r)
                               /(r*r + (self.a_s + (z*z + self.b_s*self.b_s).sqrt()).powi(2)).powf(5.0/2.0));
        let rho_thin_disk    = rho_thin_disk_z + rho_thin_disk_r;

        let rho_thick_disk_z = -self.m_g/4.0/PI*((4.0*self.a_g*z.powi(4) + (z*z + self.b_g*self.b_g).sqrt()
                               *((2.0*self.a_g*self.a_g - r*r)*z*z - self.b_g*self.b_g*r*r - 3.0*self.a_g*self.a_g*self.b_g*self.b_g)
                               + (z*z + self.b_g*self.b_g).powf(3.0/2.0)*(2.0*z*z - self.b_g*self.b_g) + self.a_g*self.b_g*self.b_g*z*z
                               - self.a_g*self.b_g*self.b_g*r*r - 3.0*self.a_g*(self.b_g).powi(4)-(self.a_g).powi(3)*self.b_g*self.b_g)
                               /((z*z + self.b_g*self.b_g).powf(3.0/2.0)*(((z*z + self.b_g*self.b_g).sqrt() + self.a_g).powi(2) + r*r).powf(5.0/2.0)));
        let rho_thick_disk_r = self.m_g/4.0/PI*((2.0*(r*r + (self.a_g + (z*z + self.b_g*self.b_g).sqrt()).powi(2)) - 3.0*r*r)
                               /(r*r + (self.a_g + (z*z + self.b_g*self.b_g).sqrt()).powi(2)).powf(5.0/2.0));
        let rho_thick_disk   = rho_thick_disk_z + rho_thick_disk_r;

        let rho_halo         = self.v_h*self.v_h*(3.0*self.a_h*self.a_h + r*r + z*z)
                               /(4.0*PI*self.g*((r*r + z*z + self.a_h*self.a_h).powi(2)));

        rho_bulge + rho_thin_disk + rho_thick_disk + rho_halo

    }

    // The z-component of the gravitational field obtained via the negative gradient of the above potential
    pub fn g_field_z(&self, r: f64, z: f64) -> f64 {

        let gfz_bulge      = -self.g*self.m_b*z/(r*r + z*z + self.a_b*self.a_b).powf(3.0/2.0);
        let gfz_thin_disk  = -self.g*self.m_b*z*(self.a_s + (z*z + self.b_s*self.b_s).sqrt())
                             / (r*r + (self.a_s + (z*z + self.b_s*self.b_s).sqrt()).powi(2)).powf(3.0/2.0)
                             / (z*z + self.b_s*self.b_s).sqrt();
        let gfz_thick_disk = -self.g*self.m_b*z*(self.a_g + (z*z + self.b_g*self.b_g).sqrt())
                             / (r*r + (self.a_g + (z*z + self.b_g*self.b_g).sqrt()).powi(2)).powf(3.0/2.0)
                             / (z*z + self.b_g*self.b_g).sqrt();
        let gfz_halo       = -self.v_h*self.v_h*z/(r*r + z*z + self.a_h*self.a_h);
        
        gfz_bulge + gfz_thin_disk + gfz_thick_disk + gfz_halo

    }

    // RK4 algorithm for the purpose of computing pressure
    pub fn rk4(&self, r: f64, z: f64, dz: f64, p: f64) -> f64 {

        let k1 = self.density(r, z)*self.g_field_z(r, z);
        let k2 = self.density(r, z + 0.5*dz)*self.g_field_z(r, z + 0.5*dz);
        let k3 = self.density(r, z + 0.5*dz)*self.g_field_z(r, z + 0.5*dz);
        let k4 = self.density(r, z + dz)*self.g_field_z(r, z + dz);

        p + (1.0/6.0)*dz*(k1 + 2.0*k2 + 2.0*k3 + k4)
}

    // Pressure from HSE condition assuming variability only in the z direction
    pub fn pressure_z(&self, r: f64, z: f64, dz: f64, p_0: f64) -> f64 {
        let mut p  = p_0;
        let mut ct = 0.0;
        while ct < z {
            let pre = self.rk4(r, ct, dz, p);
            ct += dz;
            p = pre;
        }
        p
    }

}