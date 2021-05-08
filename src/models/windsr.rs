#![allow(clippy::many_single_char_names)]
#![allow(clippy::suspicious_operation_groupings)]
//! This module solves for a time-steady hydrodynamic wind profile, including
//! a free-expansion zone initiated at a supersonic inlet, a stationary shock
//! wave (standing reverse shock) at a specified radius, and subsequent
//! subsonic decelerating flow to an outer boundary radius.
//!
//! Authors: Sagar Adhikari, Jonathan Zrake

use clap::{AppSettings, Clap};
use std::f64::consts::PI;
use std::fmt;

/// The adiabatic index
const GAMMA_LAW_INDEX: f64 = 4.0 / 3.0;

/// Speed of light in cm / s
const SPEED_OF_LIGHT: f64 = 2.99e10;

/// Numerical integration step
const DELTA_LOG_RADIUS: f64 = 1e-4;

#[derive(Clap, Debug)]
#[clap(version = "0.1.0", author = "S. Adhikari and J. Zrake")]
#[clap(setting = AppSettings::ColoredHelp)]
struct Opts {
    #[clap(long, about = "inner boundary (cm)", default_value = "1e8")]
    pub r_inner: f64,

    #[clap(long, about = "outer boundary (cm)", default_value = "1e11")]
    pub r_outer: f64,

    #[clap(long, about = "shock position (cm)", default_value = "1e10")]
    pub r_shock: f64,

    #[clap(
        short = 'l',
        long,
        about = "wind kinetic luminosity (erg / s)",
        default_value = "1e33"
    )]
    pub luminosity: f64,

    #[clap(
        short = 't',
        long,
        about = "wind terminal Lorentz factor",
        default_value = "100.0"
    )]
    pub gamma_terminal: f64,

    #[clap(
        short = 'i',
        long,
        about = "wind injection Lorentz factor",
        default_value = "10.0"
    )]
    pub gamma_launch: f64,
}

/// Holds the primitive hydrodynamic variables for a spherically symmetric
/// relativistic wind.
#[derive(Clone, Copy, Debug)]
struct Primitive {
    /// Radius (cm)
    r: f64,
    /// Radial four-velocity (gamma-beta; dimensionless)
    u: f64,
    /// Mass density (g / cm^3)
    d: f64,
    /// Specific enthalpy (cm^2 / s^2)
    h: f64,
}

impl fmt::Display for Primitive {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:.12e} {:.12e} {:.12e} {:.12e}",
            self.r, self.u, self.d, self.h
        )
    }
}

impl Primitive {
    /// Create hydrodynamic variables given the radius r (cm), gamma-beta u,
    /// mass flux (over c; g) f = rho u, and specific energy flux (over c;
    /// cm^2 / s^2), l.
    #[allow(clippy::many_single_char_names)]
    fn from_rufl(r: f64, u: f64, f: f64, l: f64) -> Self {
        let d = f / u;
        let h = l / (1.0 + u * u).sqrt();
        Self { r, u, d, h }
    }

    /// Create hydrodynamic variables given the radius r (cm), gamma-beta u,
    /// mass loss rate per steradian mdot (g / s / Sr), and wind luminosity
    /// per steradian edot (erg / s / Sr).
    #[allow(clippy::many_single_char_names)]
    fn from_ru_mdot_edot(r: f64, u: f64, mdot: f64, edot: f64) -> Self {
        let c = SPEED_OF_LIGHT;
        let d = mdot / r / r / u / c;
        let h = edot / mdot / (1.0 + u * u).sqrt();
        Self { r, u, d, h }
    }

    /// Return the wind luminosity per steradian (erg / s / Sr).
    fn luminosity(&self) -> f64 {
        self.r.powi(2) * self.l() * self.f() * SPEED_OF_LIGHT
    }

    /// Return the wind mass loss rate per steradian (g / s / Sr).
    fn mass_loss_rate(&self) -> f64 {
        self.r.powi(2) * self.f() * SPEED_OF_LIGHT
    }

    /// Return the mass flux (over c; g).
    fn f(&self) -> f64 {
        self.d * self.u
    }

    /// Return the specific momentum flux (cm^2 / s^2).
    fn k(&self) -> f64 {
        let c = SPEED_OF_LIGHT;
        let h = self.h;
        let u = self.u;
        let d = self.d;
        let mu = h - c * c;
        let e = mu / GAMMA_LAW_INDEX;
        let p = d * e * (GAMMA_LAW_INDEX - 1.0);
        h * u + p / d / u
    }

    /// Return the specific luminosity (over c; cm^2 / s^2).
    fn l(&self) -> f64 {
        let h = self.h;
        let u = self.u;
        h * (1.0 + u * u).sqrt()
    }
}

/// Return the spatial derivative du/dr for an energy and momentum-conserving
/// steady-state wind.
fn du_dr(primitive: Primitive) -> f64 {
    let Primitive { r, u, d: _, h } = primitive;
    let c = SPEED_OF_LIGHT;
    let gm = GAMMA_LAW_INDEX;
    let qh = (gm - 1.0) / gm;
    let hg = h;
    let mu = hg - c * c;
    let qd = (hg / c / c - 1.0) * qh;
    let beta = u / (1.0 + u * u).sqrt();
    let beta_f = (c * c / hg * qd / (1.0 - qh)).sqrt();
    -u / r * 2.0 * mu / hg / (qh - 1.0) * qh / (beta.powi(2) - beta_f.powi(2))
}

/// Return the jump condition for a wind with specific momentum flux (k) and
/// specific luminosity (l), and a guess value for the specific internal
/// energy (e). This function returns zero when e is either of the two
/// possible values of the specific internal energy for given k and l values.
fn jump_condition(e: f64, k: f64, l: f64) -> f64 {
    let gm = GAMMA_LAW_INDEX;
    let c = SPEED_OF_LIGHT;
    let h = c * c + gm * e;
    let lhs = l * l - (c * c + e) * h;
    let rhs = k * (l * l - h * h).sqrt();
    (lhs - rhs) / c.powi(4)
}

/// Solves for the primitive hydrodynamic state on the downstream side of a
/// standing shock wave. The input primitive state is assumed to be
/// supersonic, and the output state is subsonic.
fn solve_jump_condition(primitive: Primitive) -> Primitive {
    let gm = GAMMA_LAW_INDEX;
    let c = SPEED_OF_LIGHT;
    let f = primitive.f();
    let k = primitive.k();
    let l = primitive.l();

    let mut e1 = (l - c * c) / gm;
    let mut e0 = e1 * 0.99999;
    let mut f0 = jump_condition(e0, k, l);
    let mut f1 = jump_condition(e1, k, l);

    while f1.abs() > 1e-10 {
        let e2 = (e0 * f1 - e1 * f0) / (f1 - f0);
        let f2 = jump_condition(e2, k, l);

        e0 = e1;
        e1 = e2;
        f0 = f1;
        f1 = f2;
    }
    let e = e1;
    let h = c * c + gm * e;
    let g = l / h;
    let u = (g * g - 1.0).sqrt();
    Primitive::from_rufl(primitive.r, u, f, l)
}

fn wind_inlet(opts: &Opts) -> Primitive {
    let c = SPEED_OF_LIGHT;
    let edot = opts.luminosity / 4.0 / PI;
    let r = opts.r_inner;
    let g = opts.gamma_launch;
    let u = (g * g - 1.0).sqrt();
    let h = opts.gamma_terminal / g * c * c;
    let mdot = edot / g / h;
    let d = mdot / r / r / u / c;
    Primitive { r, u, d, h }
}

fn write_ascii_stdout<T: fmt::Display>(solution: &[T]) {
    for value in solution {
        println!("{}", value)
    }
}

fn main() {
    let opts = Opts::parse();

    let inlet_prim = wind_inlet(&opts);
    let mut r = opts.r_inner;
    let mut u = opts.gamma_launch;
    let rshock = opts.r_shock;
    let rmax = opts.r_outer;
    let edot = inlet_prim.luminosity();
    let mdot = inlet_prim.mass_loss_rate();

    let mut solution = Vec::new();

    while r < rshock {
        let dr = DELTA_LOG_RADIUS * r;
        let p = Primitive::from_ru_mdot_edot(r, u, mdot, edot);
        r += dr;
        u += du_dr(p) * dr;
        solution.push(p);
    }
    let pre_shock_prim = Primitive::from_ru_mdot_edot(r, u, mdot, edot);
    let pos_shock_prim = solve_jump_condition(pre_shock_prim);
    u = pos_shock_prim.u;

    while r < rmax {
        let p = Primitive::from_ru_mdot_edot(r, u, mdot, edot);
        let dr = DELTA_LOG_RADIUS * r;
        r += dr;
        u += du_dr(p) * dr;
        solution.push(p);
    }
    write_ascii_stdout(&solution);
}
