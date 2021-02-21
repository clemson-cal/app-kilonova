use std::f64::consts::PI;
use serde::{Serialize, Deserialize};
use crate::physics::{AnyPrimitive, LIGHT_SPEED};
use crate::traits::InitialModel;

static NOMINAL_LAUNCH_RADIUS: f64 = 1e8;
static UNIFORM_TEMPERATURE: f64 = 1e-3;




/**
 * Jet propagating through a kilonova debris cloud and surrounding relativistic
 * envelop
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JetInCloud {

    /// Mass of the merger ejecta cloud
    pub cloud_mass: f64,

    /// Time following the cloud onset when the jet begins
    pub engine_delay: f64,

    /// Duration of the engine
    pub engine_duration: f64,

    /// E / M c^2: M = cloud mass, E = isotropic-equivalent jet energy
    pub engine_strength: f64,

    /// Engine opening angle
    pub engine_theta: f64,

    /// Engine four-velocity
    pub engine_u: f64,

    /// Mass coordinate of the u=1 shell
    pub envelop_m1: f64,

    /// Beta (v/c) of the slowest envelop shell
    pub envelop_slowest_beta: f64,

    /// Beta (v/c) of the outer shell
    pub envelop_fastest_beta: f64,

    /// Index psi in u(m) ~ m^-psi
    pub envelop_psi: f64,
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
        self.print(&mut std::io::stdout());
        Ok(())
    }

    fn primitive_at(&self, coordinate: (f64, f64), t: f64) -> AnyPrimitive {
        let (r, q) = coordinate;
        let f = self.mass_rate_per_steradian(r, q, t);
        let u = self.gamma_beta(r, q, t);
        let d = f / (r * r * u) / LIGHT_SPEED;
        let p = d * UNIFORM_TEMPERATURE;

        AnyPrimitive{
            velocity_r: u,
            velocity_q: 0.0,
            mass_density: d,
            gas_pressure: p,
        }
    }

    fn scalar_at(&self, coordinate: (f64, f64), t: f64) -> f64 {
        let (r, q) = coordinate;

        match self.zone(r, q, t) {
            Zone::Cloud       => 1e+0,
            Zone::Jet         => 1e+2,
            Zone::Envelop     => 1e-2,
        }
    }
}




// ============================================================================
impl JetInCloud
{

    /**
     * Print a summary of the relevant times in the model.
     */
    pub fn print<Writer: std::io::Write>(&self, writer: &mut Writer) {
        writeln!(writer,
            r#"
        jet_in_cloud model description:
            t1 (slowest envelop shell comes to r=10^8 cm) = {:.04}
            t2 (jet turns on)                             = {:.04}
            t3 (jet head comes through r=10^8 cm)         = {:.04}
            t4 (jet turns off)                            = {:.04}
        "#, self.get_t1(),
            self.get_t2(),
            self.get_t3(),
            self.get_t4()).unwrap();
    }

    /**
     * The time when the slowest envelop shell comes through the launch radius
     */
    pub fn get_t1(&self) -> f64 {
        NOMINAL_LAUNCH_RADIUS / self.envelop_slowest_beta / LIGHT_SPEED
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
        self.get_t2() + NOMINAL_LAUNCH_RADIUS / self.engine_beta() / LIGHT_SPEED
    }

    /**
     * Time when the jet turns off
     */
    pub fn get_t4(&self) -> f64 {
        self.get_t2() + self.engine_duration
    }

    /**
     * Four-velocity gamma-beta of the slowest envelop shell
     */
    pub fn envelop_slowest_u(&self) -> f64 {
        let b = self.envelop_slowest_beta;
        b / (1.0 - b * b).sqrt()
    }

    /**
     * Dimensionless jet velocity: v_jet / c
     */
    pub fn engine_beta(&self) -> f64 {
        self.engine_u / (1.0 + self.engine_u.powi(2)).sqrt()
    }

    /**
     * Determine if a polar angle is within theta_jet of either pole.
     *
     * * `q` - The polar angle theta
     */
    pub fn in_nozzle(&self, q: f64) -> bool {
        q < self.engine_theta || q > PI - self.engine_theta
    }

    /**
     * Determine the zone of the ambient medium for a given radius and time.
     *
     * * `r` - Radius
     * * `q` - Polar angle
     * * `t` - Time
     */
    pub fn zone(&self, r: f64, q: f64, t: f64) -> Zone {
        let v_min = self.envelop_slowest_beta * LIGHT_SPEED;
        let v_jet = self.engine_beta() * LIGHT_SPEED;

        let r_cloud_envelop_interface = v_min * t;
        let r_jet_head = v_jet * (t - self.engine_delay);
        let r_jet_tail = v_jet * (t - self.engine_delay - self.engine_duration);

        if self.in_nozzle(q) && r < r_jet_head  && r > r_jet_tail {
            Zone::Jet
        } else if r > r_cloud_envelop_interface {
            Zone::Envelop
        } else {
            Zone::Cloud
        }
    }

    /**
     * Return the radial four-velocity (gamma-beta).
     *
     * * `r` - The radius
     * * `q` - The polar angle theta
     * * `t` - The time
     */
    pub fn gamma_beta(&self, r: f64, q: f64, t: f64) -> f64 {
        match self.zone(r, q, t) {
            Zone::Cloud => {
                self.envelop_slowest_u()
            }
            Zone::Envelop => {
                let b = f64::min(r / t / LIGHT_SPEED, self.envelop_fastest_beta);
                let u = b / f64::sqrt(1.0 - b * b);
                u
            }
            Zone::Jet => {
                self.engine_u
            }
        }
    }

    /**
     * Return the mass rate per solid angle.
     *
     * * `r` - The radius
     * * `q` - The polar angle theta
     * * `t` - The time
     */
    pub fn mass_rate_per_steradian(&self, r: f64, q: f64, t: f64) -> f64 {

        match self.zone(r, q, t) {
            Zone::Cloud => {
                self.cloud_mass_rate_per_steradian()
            }
            Zone::Envelop => {
                let s = f64::min(r / t / LIGHT_SPEED, self.envelop_fastest_beta);
                let f = f64::powf(s, -1.0 / self.envelop_psi) * f64::powf(1.0 - s * s, 0.5 / self.envelop_psi - 1.0);
                self.envelop_m1 / (4.0 * PI * self.envelop_psi * t) * f
            }
            Zone::Jet => {
                self.jet_mass_rate_per_steradian()
            }
        }
    }

    // fn taper(x: f64)
    fn jet_mass_rate_per_steradian(&self) -> f64 {
        let engine_gamma = f64::sqrt(1.0 + self.engine_u * self.engine_u);
        let e = self.engine_strength * self.cloud_mass;
        let l = e / (4.0 * PI * self.engine_duration);
        l / engine_gamma
    }

    fn cloud_mass_rate_per_steradian(&self) -> f64 {
        self.cloud_mass / (4.0 * PI * self.engine_delay)
    }
}
