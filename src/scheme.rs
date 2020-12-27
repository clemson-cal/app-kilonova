use crate::mesh::Mesh;
use crate::state::State;
use crate::traits::{Conserved, Hydrodynamics};




// ============================================================================
#[allow(unused)]
pub fn advance<H, C>(state: &mut State<C>, hydro: &H, mesh: &Mesh)
where
	H: Hydrodynamics<Conserved = C>,
	C: Conserved {
	state.time += 0.1;
	state.iteration += 1;
	std::thread::sleep(std::time::Duration::from_millis(100));
}
