use std::f64::consts::PI;
use serde::{Serialize, Deserialize};
use crate::physics::{AgnosticPrimitive, LIGHT_SPEED};
use crate::traits::InitialModel;

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




/**
 * Different space-time zones in the setup
 */
pub enum Zone {
    Envelop,
    Cloud,
    Jet,
}




// ============================================================================
impl InitialModel for JetInCloud {

    fn validate(&self) -> anyhow::Result<()> {
        self.print();
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

        match self.zone(r, q, t) {
            Zone::Cloud => self.envelop_mass * 1e2,
            Zone::Envelop => self.envelop_mass * (self.gamma_beta(r, q, t) / self.u_min).powf(-1.0 / self.psi),
            Zone::Jet => self.envelop_mass * 1e-12
        }
    }
}




// ============================================================================
impl JetInCloud
{

    pub fn print(&self) {
        println!("\tjet_in_cloud model description:\n");
        println!("\tt1 (slowest envelop shell comes through launch_radius) = {:.04}", self.get_t1());
        println!("\tt2 (the jet turns on)                                  = {:.04}", self.get_t2());
        println!("\tt3 (jet head comes through launch_radius)              = {:.04}", self.get_t3());
        println!("\tt4 (time when the jet turns off)                       = {:.04}", self.get_t4());
        println!();
    }

    pub fn engine_beta(&self) -> f64 {
        self.engine_u / (1.0 + self.engine_u.powi(2)).sqrt()
    }

    /**
     * The time when the slowest envelop shell comes through the launch radius
     */
    pub fn get_t1(&self) -> f64 {
        let v_min = self.u_min / (1.0 + self.u_min.powi(2)).sqrt() * LIGHT_SPEED;
        self.launch_radius / v_min
    }

    /**
     * Time when the jet turns on
     */
    pub fn get_t2(&self) -> f64 {
        self.engine_delay
    }

    /**
     * Time when the jet comes through the lauch radius
     */
    pub fn get_t3(&self) -> f64 {
        self.get_t2() + self.launch_radius / self.engine_beta() / LIGHT_SPEED
    }

    /**
     * Time when the jet turns off
     */
    pub fn get_t4(&self) -> f64 {
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
     * Determine the zone of the ambient medium for a given radius and time:
     *
     * * `r` - Radius
     * * `q` - Polar angle
     * * `t` - Time
     */
    pub fn zone(&self, r: f64, q: f64, t: f64) -> Zone {
        let v_min = self.u_min / f64::sqrt(1.0 + self.u_min * self.u_min) * LIGHT_SPEED;
        let v_jet = self.engine_beta() * LIGHT_SPEED;

        let r_cloud_envelop_interface = v_min * t;
        let r_jet_head = v_jet * (t - self.engine_delay);
        let r_jet_tail = v_jet * (t - self.engine_delay - self.engine_duration);

        if self.in_nozzle(q) && r < r_jet_head && r > r_jet_tail {
            Zone::Jet
        } else if r > r_cloud_envelop_interface {
            Zone::Envelop
        } else {
            Zone::Cloud
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

        match self.zone(r, q, t) {
            Zone::Cloud => {
                self.u_min
            },
            Zone::Envelop => {
                let b = f64::min(r / t / LIGHT_SPEED, MAX_BETA);
                let u = b / f64::sqrt(1.0 - b * b);
                u
            },
            Zone::Jet => {
                self.engine_u
            }
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

        match self.zone(r, q, t) {
            Zone::Cloud => {
                mc / (4.0 * PI * self.engine_delay)
            },
            Zone::Envelop => {
                let s = f64::min(r / t / LIGHT_SPEED, MAX_BETA);
                let f = f64::powf(s / self.u_min, -1.0 / self.psi) * f64::powf(1.0 - s * s, 0.5 / self.psi - 1.0);
                m0 / (4.0 * PI * self.psi * t) * f
            },
            Zone::Jet => {
                self.get_engine_mass_flux()                
            }
        }
    }
}
