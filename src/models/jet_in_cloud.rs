use std::f64::consts::PI;
use serde::{Serialize, Deserialize};
use crate::physics::AgnosticPrimitive;
use crate::traits::InitialModel;

static LIGHT_SPEED: f64 = 3e10;
static UNIFORM_ENTROPY: f64 = 1e-4;
static MAX_BETA: f64 = 0.97;
static GAMMA_LAW_INDEX: f64 = 4.0 / 3.0;




/**
 * Jet propagating through a kilonova debris cloud and surrounding relativistic
 * envelop
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JetInCloud {

    /// Four-velocity of the slowest envelop shell
    pub launch_radius: f64,

    /// Index psi in u(m) ~ m^-psi
    pub cloud_mass: f64,

    /// Mass of the relativistic envelop
    pub engine_delay: f64,

    /// Engine four-velocity
    pub engine_duration: f64,

    /// Engine opening angle
    pub engine_strength: f64,

    /// E / M c^2: M = cloud mass, E = isotropic-equivalent jet energy
    pub engine_theta: f64,

    /// Duration of the engine
    pub engine_u: f64,

    /// Time following the cloud onset when the jet begins
    pub envelop_mass: f64,

    /// Mass of the merger ejecta cloud
    pub psi: f64,

    /// Radius where the inflow starts from [cm]
    pub u_min: f64,
}




// ============================================================================
impl InitialModel for JetInCloud {

    fn validate(&self) -> anyhow::Result<()> {
        Ok(())
    }

    fn primitive_at(&self, coordinate: (f64, f64), t: f64) -> AgnosticPrimitive {
        let (r, q) = coordinate;
        let f = self.mass_flux(r, q, t);
        let u = self.gamma_beta(r, q, t);
        let d = f / (r * r * u);
        let d0 = self.cloud_mass / self.launch_radius.powi(3);
        let s = UNIFORM_ENTROPY;
        let p = s * f64::powf(d / d0, GAMMA_LAW_INDEX); // TODO: load gamma from hydro

        AgnosticPrimitive{
            velocity_r: u,
            velocity_q: 0.0,
            mass_density: d,
            gas_pressure: p,
        }
    }

    fn scalar_at(&self, coordinate: (f64, f64), t: f64) -> f64 {
        let (r, q) = coordinate;

        if self.get_zone(r, t) == 1 {
            let u = self.gamma_beta(r, q, t);
            self.envelop_mass * (u / self.u_min).powf(-1.0 / self.psi)
        } else if self.get_zone(r, t) == 2 {
            self.envelop_mass * 0.5
        } else if self.get_zone(r, t) == 3 && self.in_nozzle(q) {
            1e-12
        } else {
            self.envelop_mass * 0.5
        }
    }
}




// ============================================================================
impl JetInCloud
{

    pub fn print(&self) {
        println!("\tjet-cloud model description:\n\n");
        println!("\tt1 (time when slowest envelop shell comes through r=1) = {}", self.get_t1());
        println!("\tt2 (time when the jet turns on)                        = {}", self.get_t2());
        println!("\tt3 (time when the jet turns off)                       = {}", self.get_t3());
        println!();
    }


    /**
     * The time when the slowest envelop shell comes through the launch radius
     */
    pub fn get_t1(&self) -> f64 {
        let v_min = self.u_min / f64::sqrt(1.0 + self.u_min * self.u_min) * LIGHT_SPEED;
        self.launch_radius / v_min
    }


    /**
     * The time when the jet turns on
     */
    pub fn get_t2(&self) -> f64 {
        self.get_t1() + self.engine_delay
    }


    /**
     * The time when the jet turns off
     */
    pub fn get_t3(&self) -> f64 {
        self.get_t2() + self.engine_duration
    }


    /**
     * Return the mass flux, per sterian, in the jet
     */
    pub fn get_engine_mass_flux(&self) -> f64 {
        let engine_gamma = f64::sqrt(1.0 + self.engine_u * self.engine_u);
        let e = self.engine_strength * self.cloud_mass;
        let l = e / (4.0 * PI * self.engine_duration);
        let mdot = l / engine_gamma;
        mdot
    }


    /**
     * Return the true of a polar angle is within theta_jet of either pole
     *
     * * `q` - The polar angle theta
     */
    pub fn in_nozzle(&self, q: f64) -> bool {
        q < self.engine_theta || q > PI - self.engine_theta
    }


    /**
     * Determine the zone of the ambient medium for a given radius and time
     *
     * * `r` - The radius
     * * `t` - The time
     */
    pub fn get_zone(&self, r: f64, t: f64) -> usize {
        let v_min = self.u_min / f64::sqrt(1.0 + self.u_min * self.u_min) * LIGHT_SPEED;

        if t <= r / v_min {
            1
        }
        else if t <= r / v_min + self.engine_delay {
            2
        }
        else if t <= r / v_min + self.engine_delay + self.engine_duration {
            3
        }
        else {
            4
        }
    }


    /**
     * Return the radial four-velocity (gamma-beta)
     *
     * * `r` - The radius
     * * `q` - The polar angle theta
     * * `t` - The time
     */
    pub fn gamma_beta(&self, r: f64, q: f64, t: f64) -> f64 {
        if self.get_zone(r, t) == 1 {
            let b = f64::min(r / t / LIGHT_SPEED, MAX_BETA);
            let u = b / f64::sqrt(1.0 - b * b);
            u
        }
        else if self.get_zone(r, t) == 2 {
            self.u_min
        }
        else if (self.get_zone(r, t) == 3 || self.get_zone(r, t) == 4) && self.in_nozzle(q) {
            if t > self.get_t2() && t < self.get_t3() {
                self.engine_u
            } else {
                self.u_min
            }
        }
        else {
            self.u_min
        }
    }


    /**
     * Return the mass flux per solid angle
     *
     * * `r` - The radius
     * * `q` - The polar angle theta
     * * `t` - The time
     */
    pub fn mass_flux(&self, r: f64, q: f64, t: f64) -> f64 {
        let m0 = self.envelop_mass;
        let mc = self.cloud_mass;

        if self.get_zone(r, t) == 1 {
            let s = f64::min(r / t / LIGHT_SPEED, MAX_BETA);
            let f = f64::powf(s / self.u_min, -1.0 / self.psi) * f64::powf(1.0 - s * s, 0.5 / self.psi - 1.0);
            m0 / (4.0 * PI * self.psi * t) * f
        }
        else if self.get_zone(r, t) == 2 {
            mc / (4.0 * PI * self.engine_delay)
        }
        else if self.get_zone(r, t) == 3 && self.in_nozzle(q) {
            self.get_engine_mass_flux()
        }
        else {
            mc / (4.0 * PI * self.engine_delay)
        }
    }
}
