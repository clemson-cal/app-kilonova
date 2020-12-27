use serde::{Serialize, Deserialize};
use crate::physics::AgnosticPrimitive;
use crate::traits::InitialModel;




/**
 * Jet propagating through a kilonova debris cloud and surrounding relativistic
 * envelop
 */
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JetInCloud {
    pub engine_duration: f64,
    pub engine_strength: f64,
}




// ============================================================================
impl InitialModel for JetInCloud {
    fn primitive_at(&self, _coordinate: (f64, f64)) -> AgnosticPrimitive {
        AgnosticPrimitive{
            velocity_r: 0.0,
            velocity_q: 0.0,
            mass_density: 1.0,
            gas_pressure: 0.01,
        }
    }

    fn scalar_at(&self, _coordinate: (f64, f64)) -> f64 {
        0.0
    }
}
