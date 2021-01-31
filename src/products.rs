use std::collections::HashMap;
use serde::{Serialize, Deserialize};
use ndarray::{ArcArray, Ix1, Ix2};
use crate::app::{self, Configuration, AgnosticHydro, AgnosticState};
use crate::mesh::{BlockIndex, GridGeometry};
use crate::physics::AgnosticPrimitive;
use crate::products;
use crate::state::{BlockState, State};
use crate::traits::{Conserved, Hydrodynamics};




/**
 * Useful per-block data for post-processing and plotting
 */
#[derive(Clone, Serialize, Deserialize)]
pub struct BlockProducts {
	pub radial_vertices: ArcArray<f64, Ix1>,
	pub polar_vertices: ArcArray<f64, Ix1>,
	pub primitive: ArcArray<AgnosticPrimitive, Ix2>,
	pub scalar: ArcArray<f64, Ix2>,	
}




/**
 * Useful data for post-processing and plotting
 */
#[derive(Serialize, Deserialize)]
pub struct Products {
	pub time: f64,
	pub blocks: HashMap<BlockIndex, BlockProducts>,
	pub config: Configuration,
	pub version: String,
}




// ============================================================================
impl BlockProducts {
	pub fn from_block_state<H, C>(state: &BlockState<C>, hydro: &H, geometry: &GridGeometry) -> Self
	where
		H: Hydrodynamics<Conserved = C>,
		C: Conserved {

		let scalar = &state.scalar_mass / &state.conserved.mapv(|u| u.lab_frame_mass());
		let primitive = (&state.conserved / &geometry.cell_volumes)
			.mapv(|q| hydro.to_primitive(q))
			.mapv(|p| hydro.agnostic(&p));

		BlockProducts{
			radial_vertices: geometry.radial_vertices.clone(),
			polar_vertices: geometry.polar_vertices.clone(),
			primitive: primitive.to_shared(),
			scalar: scalar.to_shared(),
		}
	}
}




// ============================================================================
impl Products {
	pub fn from_state<H, C>(state: &State<C>, hydro: &H, config: &Configuration) -> Self
	where
		H: Hydrodynamics<Conserved = C>,
		C: Conserved {

		let geometry = config.mesh.grid_blocks_geometry(state.time);
		let mut blocks = HashMap::new();

		for (index, block_state) in &state.solution {
			blocks.insert(*index, BlockProducts::from_block_state(block_state, hydro, &geometry[index]));
		}

		Products{
			time: state.time,
			blocks: blocks,
			config: config.clone(),
			version: app::VERSION_AND_BUILD.to_string(),
		}
	}
	pub fn from_app(app: &app::App) -> Self {
		match (&app.state, &app.config.hydro) {
			(AgnosticState::Newtonian(state), AgnosticHydro::Newtonian(hydro)) => {
				products::Products::from_state(state, hydro, &app.config)
			},
			(AgnosticState::Relativistic(state), AgnosticHydro::Relativistic(hydro)) => {
				products::Products::from_state(state, hydro, &app.config)
			},
			_ => unreachable!()
		}
	}
}
