use std::collections::HashMap;
use serde::{Serialize, Deserialize};
use ndarray::{ArcArray, Ix1, Ix2};
use crate::Configuration;
use crate::traits::{Conserved, Hydrodynamics, Primitive};
use crate::state::{BlockState, State};
use crate::mesh::{BlockIndex, GridGeometry, Mesh};




/**
 * Useful per-block data for post-processing and plotting
 */
#[derive(Serialize, Deserialize)]
pub struct BlockProducts<P: Primitive> {
	pub radial_vertices: ArcArray<f64, Ix1>,
	pub polar_vertices: ArcArray<f64, Ix1>,
	pub primitive: ArcArray<P, Ix2>,
	pub scalar: ArcArray<f64, Ix2>,	
}




/**
 * Useful data for post-processing and plotting
 */
#[derive(Serialize, Deserialize)]
pub struct Products<P: Primitive> {
	time: f64,
	blocks: HashMap<BlockIndex, BlockProducts<P>>,
	config: Configuration,
	version: String,
}




// ============================================================================
impl<P: Primitive> BlockProducts<P> {
	pub fn from_block_state<H, C>(state: &BlockState<C>, hydro: &H, geometry: &GridGeometry) -> anyhow::Result<Self>
	where
		H: Hydrodynamics<Conserved = C, Primitive = P>,
		C: Conserved {

		let scalar = &state.scalar_mass / &state.conserved.mapv(|u| u.lab_frame_mass());
		let primitive = (&state.conserved / &geometry.cell_volumes).mapv(|q| hydro.to_primitive(q));

		Ok(BlockProducts{
			radial_vertices: geometry.radial_vertices.clone(),
			polar_vertices: geometry.polar_vertices.clone(),
			primitive: primitive.to_shared(),
			scalar: scalar.to_shared(),
		})
	}
}




// ============================================================================
impl<P: Primitive> Products<P> {
	pub fn from_state<H, C>(state: &State<C>, hydro: &H, mesh: &Mesh, config: &Configuration) -> anyhow::Result<Self>
	where
		H: Hydrodynamics<Conserved = C, Primitive = P>,
		C: Conserved {

		let geometry = mesh.grid_blocks_geometry();
		let mut blocks = HashMap::new();

		for (index, block_state) in &state.solution {
			blocks.insert(*index, BlockProducts::from_block_state(block_state, hydro, &geometry[index])?);
		}

		Ok(Products{
			time: state.time,
			blocks: blocks,
			config: config.clone(),
			version: crate::VERSION_AND_BUILD.to_string(),
		})
	}
}
