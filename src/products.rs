use std::collections::HashMap;
use serde::{Serialize, Deserialize};
use ndarray::{ArcArray, Ix1, Ix2};
use crate::scheme::{try_block_primitive};
use crate::app::{self, Configuration, AnyHydro, AnyState};
use crate::mesh::{BlockIndex, GridGeometry};
use crate::physics::{AnyPrimitive};
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
	pub primitive: ArcArray<AnyPrimitive, Ix2>,
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
		let primitive= {
			try_block_primitive(hydro, state.conserved.to_shared() / &geometry.cell_volumes, &geometry)
				.map(|p| p.to_shared())
				.unwrap()
				.mapv(|p| hydro.any(&p))
		};

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
			(AnyState::Newtonian(state), AnyHydro::Newtonian(hydro)) => {
				products::Products::from_state(state, hydro, &app.config)
			},
			(AnyState::Relativistic(state), AnyHydro::Relativistic(hydro)) => {
				products::Products::from_state(state, hydro, &app.config)
			},
			_ => unreachable!()
		}
	}
}
