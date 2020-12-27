use serde::{Serialize, Deserialize};
use crate::traits::InitialModel;
use crate::physics::AgnosticPrimitive;




/**
 * Explosion in a horizontally stratified external medium
 */
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HaloKilonova {
    pub explosion_energy: f64,
}




// ============================================================================
impl InitialModel for HaloKilonova {
    fn at(&self, _coordinate: (f64, f64)) -> AgnosticPrimitive {
    	AgnosticPrimitive{
    		velocity_r: 0.0,
    		velocity_q: 0.0,
    		mass_density: 1.0,
    		gas_pressure: 0.01,
    	}
    }
}
