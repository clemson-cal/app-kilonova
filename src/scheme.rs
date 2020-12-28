#![allow(unused)]
use std::collections::HashMap;
use ndarray::{ArcArray, Ix2, Axis, stack, s};
use crate::Model;
use crate::mesh::{GridGeometry, Mesh};
use crate::state::{State, BlockState};
use crate::traits::{Conserved, Hydrodynamics, Primitive};




// ============================================================================
#[allow(unused)]
pub fn advance<H, C>(state: &mut State<C>, hydro: &H, model: &Model, mesh: &Mesh) -> anyhow::Result<()>
where
	H: Hydrodynamics<Conserved = C>,
	C: Conserved {

	let block_geometry = mesh.grid_blocks_geometry()?;
	let mut primitive_map = HashMap::new();
	let mut scalar_map = HashMap::new();

	let mut insert_into_map = |index, state: &BlockState<C>, geometry: &GridGeometry| {
		let scalar = &state.scalar_mass / &state.conserved.mapv(|u| u.lab_frame_mass());
		let primitive = (&state.conserved / &geometry.cell_volumes).mapv(|q| hydro.to_primitive(q));
		primitive_map.insert(index, primitive);
		scalar_map.insert(index, scalar);
	};

	for (&index, state) in &state.solution {
		insert_into_map(index, state, &block_geometry[&index]);
	}

	let (inner_bnd_index, outer_bnd_index) = state.inner_outer_boundary_indexes();

	let inner_bnd_geom = mesh.subgrid(inner_bnd_index)?.geometry();
	let outer_bnd_geom = mesh.subgrid(outer_bnd_index)?.geometry();

	let inner_bnd_state = BlockState::from_model(model, hydro, &inner_bnd_geom, state.time);
	let outer_bnd_state = BlockState::from_model(model, hydro, &outer_bnd_geom, state.time);

	insert_into_map(inner_bnd_index, &inner_bnd_state, &inner_bnd_geom);
	insert_into_map(outer_bnd_index, &outer_bnd_state, &outer_bnd_geom);

	for (index, state) in &state.solution {
		let il = (index.0 - 1, index.1);
		let i0 = (index.0,     index.1);
		let ir = (index.0 + 1, index.1);

		let pl = &primitive_map[&il];
		let p0 = &primitive_map[&i0];
		let pr = &primitive_map[&ir];
		let pe = stack(Axis(0), &[pl.slice(s![-2.., ..]), p0.view(), pr.slice(s![..2, ..])]);

		let sl = &scalar_map[&il];
		let s0 = &scalar_map[&i0];
		let sr = &scalar_map[&ir];
		let se = stack(Axis(0), &[sl.slice(s![-2.., ..]), s0.view(), sr.slice(s![..2, ..])]);
	}

	state.time += 0.1;
	state.iteration += 1;
	std::thread::sleep(std::time::Duration::from_millis(100));

	Ok(())
}
