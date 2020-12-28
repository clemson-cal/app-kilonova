use serde::{Serialize, Deserialize};
use crate::traits::InitialModel;
use crate::physics::AgnosticPrimitive;




/**
 * Explosion in a horizontally stratified external medium
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HaloKilonova {
    pub explosion_energy: f64,
}




// ============================================================================
impl InitialModel for HaloKilonova {
    fn primitive_at(&self, _coordinate: (f64, f64), _time: f64) -> AgnosticPrimitive {
    	AgnosticPrimitive{
    		velocity_r: 0.0,
    		velocity_q: 0.0,
    		mass_density: 1.0,
    		gas_pressure: 0.01,
    	}
    }

    fn scalar_at(&self, _coordinate: (f64, f64), _time: f64) -> f64 {
        0.0
    }
}
