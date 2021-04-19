//! 10-parameter axisymmetric model of a disk galaxy, using the Plummer model
//! for the central bulge, two Miyamoto- Nagai models for the thin and thick
//! disks, and a logarithmic model for the dark matter halo.
//!

use std::f64::consts::PI;

/// Galactic model parameters, including the gravitational constant. Here, slr
/// stands for solar masses.
/// 
#[derive(Clone, Debug)]
pub struct GalacticModel {
    pub g: f64,   // gravitational constant (kpc*kpc*kpc/Myr/Myr/slr)
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

/// The components of the galactic model: these can be the model contributions
/// to the mass density, gravitational potential, or gravitational
/// acceleration, depending on the context.
///
pub struct ModelComponents {
    pub bulge: f64,
    pub thin_disk: f64,
    pub thick_disk: f64,
    pub halo: f64,
}

impl ModelComponents {
    pub fn total(&self) -> f64 {
        self.bulge + self.thin_disk + self.thick_disk + self.halo
    }
}

impl GalacticModel {
    /// Gravitational potential
    /// 
    pub fn potential(&self, r: f64, z: f64) -> ModelComponents {
        let Self {
            g,
            m_b,
            a_b,
            v_h,
            a_h,
            m_s,
            a_s,
            b_s,
            m_g,
            a_g,
            b_g,
        } = self.clone();

        let bulge = -g * m_b / (r * r + z * z + a_b * a_b).sqrt();
        let thin_disk = -g * m_s / (r * r + (a_s + (z * z + b_s * b_s).sqrt()).powi(2)).sqrt();
        let thick_disk = -g * m_g / (r * r + (a_g + (z * z + b_g * b_g).sqrt()).powi(2)).sqrt();
        let halo = 0.5 * v_h * v_h * (r * r + z * z + a_h * a_h).ln();

        ModelComponents {
            bulge,
            thin_disk,
            thick_disk,
            halo,
        }
    }

    /// Mass density profile obtained via Poisson's equation using the above
    /// potential.
    /// 
    pub fn density(&self, r: f64, z: f64) -> ModelComponents {
        let Self {
            g,
            m_b,
            a_b,
            v_h,
            a_h,
            m_s,
            a_s,
            b_s,
            m_g,
            a_g,
            b_g,
        } = self.clone();

        let bulge = (3.0 * m_b / 4.0 / PI / a_b / a_b / a_b)
            * (1.0 + (r * r + z * z) / a_b / a_b).powf(-2.5);

        let thin_disk_z = -m_s / 4.0 / PI
            * ((4.0 * a_s * z.powi(4)
                + (z * z + b_s * b_s).sqrt()
                    * ((2.0 * a_s * a_s - r * r) * z * z
                        - b_s * b_s * r * r
                        - 3.0 * a_s * a_s * b_s * b_s)
                + (z * z + b_s * b_s).powf(1.5) * (2.0 * z * z - b_s * b_s)
                + a_s * b_s * b_s * z * z
                - a_s * b_s * b_s * r * r
                - 3.0 * a_s * b_s.powi(4)
                - a_s.powi(3) * b_s * b_s)
                / ((z * z + b_s * b_s).powf(1.5)
                    * (((z * z + b_s * b_s).sqrt() + a_s).powi(2) + r * r).powf(2.5)));

        let thin_disk_r = m_s / 4.0 / PI
            * ((2.0 * (r * r + (a_s + (z * z + b_s * b_s).sqrt()).powi(2)) - 3.0 * r * r)
                / (r * r + (a_s + (z * z + b_s * b_s).sqrt()).powi(2)).powf(2.5));

        let thick_disk_z = -m_g / 4.0 / PI
            * ((4.0 * a_g * z.powi(4)
                + (z * z + b_g * b_g).sqrt()
                    * ((2.0 * a_g * a_g - r * r) * z * z
                        - b_g * b_g * r * r
                        - 3.0 * a_g * a_g * b_g * b_g)
                + (z * z + b_g * b_g).powf(1.5) * (2.0 * z * z - b_g * b_g)
                + a_g * b_g * b_g * z * z
                - a_g * b_g * b_g * r * r
                - 3.0 * a_g * b_g.powi(4)
                - a_g.powi(3) * b_g * b_g)
                / ((z * z + b_g * b_g).powf(1.5)
                    * (((z * z + b_g * b_g).sqrt() + a_g).powi(2) + r * r).powf(2.5)));

        let thick_disk_r = m_g / 4.0 / PI
            * ((2.0 * (r * r + (a_g + (z * z + b_g * b_g).sqrt()).powi(2)) - 3.0 * r * r)
                / (r * r + (a_g + (z * z + b_g * b_g).sqrt()).powi(2)).powf(2.5));

        let halo = v_h * v_h * (3.0 * a_h * a_h + r * r + z * z)
            / (4.0 * PI * g * (r * r + z * z + a_h * a_h).powi(2));

        let thin_disk = thin_disk_z + thin_disk_r;
        let thick_disk = thick_disk_z + thick_disk_r;

        ModelComponents {
            bulge,
            thin_disk,
            thick_disk,
            halo,
        }
    }

    /// The z-component of the gravitational field obtained via the negative
    /// gradient of the above potential.
    /// 
    pub fn g_field_z(&self, r: f64, z: f64) -> ModelComponents {
        let Self {
            g,
            m_b,
            a_b,
            v_h,
            a_h,
            m_s: _,
            a_s,
            b_s,
            m_g: _,
            a_g,
            b_g,
        } = self.clone();

        let bulge = -g * m_b * z / (r * r + z * z + a_b * a_b).powf(1.5);

        let thin_disk = -g * m_b * z * (a_s + (z * z + b_s * b_s).sqrt())
            / (r * r + (a_s + (z * z + b_s * b_s).sqrt()).powi(2)).powf(1.5)
            / (z * z + b_s * b_s).sqrt();

        let thick_disk = -g * m_b * z * (a_g + (z * z + b_g * b_g).sqrt())
            / (r * r + (a_g + (z * z + b_g * b_g).sqrt()).powi(2)).powf(3.0 / 2.0)
            / (z * z + b_g * b_g).sqrt();

        let halo = -v_h * v_h * z / (r * r + z * z + a_h * a_h);

        ModelComponents {
            bulge,
            thin_disk,
            thick_disk,
            halo,
        }
    }

    /// RK4 algorithm for the purpose of computing pressure.
    ///
    pub fn pressure_difference_rk4(&self, r: f64, z: f64, dz: f64) -> f64 {
        let rho = |r, z| self.density(r, z).thin_disk;
        let gz = |r, z| self.g_field_z(r, z).total();

        let k1 = rho(r, z) * gz(r, z);
        let k2 = rho(r, z + 0.5 * dz) * gz(r, z + 0.5 * dz);
        let k3 = rho(r, z + 0.5 * dz) * gz(r, z + 0.5 * dz);
        let k4 = rho(r, z + dz) * gz(r, z + dz);

        dz * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    }

    /// Pressure of the gas in the thin disk from HSE condition, assuming the
    /// pressure varies only in the z direction. Returns a `Vec` where each
    /// element is a two-tuple of (altitude, pressure).
    ///
    pub fn vertical_pressure_profile(&self, r: f64, zmax: f64, dz: f64, mut p: f64) -> Vec<(f64, f64)> {
        let mut profile = Vec::new();
        let mut z = 0.0;
        while z < zmax {
            profile.push((z, p));
            p += self.pressure_difference_rk4(r, z, dz);
            z += dz;
        }
        profile.push((z, p));
        profile
    }
}