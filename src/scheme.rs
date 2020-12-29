#![allow(unused)]
use std::collections::HashMap;
use ndarray::{ArcArray, Array, Ix2, Axis, concatenate, s};
use crate::mesh::{BlockIndex, GridGeometry, Mesh};
use crate::Model;
use crate::physics::Direction;
use crate::state::{State, BlockState};
use crate::traits::{Conserved, Hydrodynamics, Primitive};

static PLM_THETA: f64 = 1.5;




// ============================================================================
#[allow(unused)]
pub fn advance<H, C>(state: &mut State<C>, hydro: &H, model: &Model, mesh: &Mesh, geometry: &HashMap<BlockIndex, GridGeometry>)
	-> anyhow::Result<()>
where
	H: Hydrodynamics<Conserved = C>,
	C: Conserved {

	let mut primitive_map = HashMap::new();
	let mut scalar_map = HashMap::new();

	let mut insert_into_map = |index, state: &BlockState<C>, geometry: &GridGeometry| {
		let scalar = &state.scalar_mass / &state.conserved.mapv(|u| u.lab_frame_mass());
		let primitive = (&state.conserved / &geometry.cell_volumes).mapv(|q| hydro.to_primitive(q));
		primitive_map.insert(index, primitive);
		scalar_map.insert(index, scalar);
	};

	for (&index, state) in &state.solution {
		insert_into_map(index, state, &geometry[&index]);
	}

	let (inner_bnd_index, outer_bnd_index) = state.inner_outer_boundary_indexes();

	let inner_bnd_geom = mesh.subgrid(inner_bnd_index).geometry();
	let outer_bnd_geom = mesh.subgrid(outer_bnd_index).geometry();

	let inner_bnd_state = BlockState::from_model(model, hydro, &inner_bnd_geom, state.time);
	let outer_bnd_state = BlockState::from_model(model, hydro, &outer_bnd_geom, state.time);

	insert_into_map(inner_bnd_index, &inner_bnd_state, &inner_bnd_geom);
	insert_into_map(outer_bnd_index, &outer_bnd_state, &outer_bnd_geom);

	for (index, state) in &state.solution {
		use godunov_core::piecewise_linear::plm_gradient;

		let il = (index.0 - 1, index.1);
		let i0 = (index.0,     index.1);
		let ir = (index.0 + 1, index.1);

		let pl = &primitive_map[&il];
		let p0 = &primitive_map[&i0];
		let pr = &primitive_map[&ir];
		let pe = concatenate(Axis(0), &[pl.slice(s![-2.., ..]), p0.view(), pr.slice(s![..2, ..])])?;

		let sl = &scalar_map[&il];
		let s0 = &scalar_map[&i0];
		let sr = &scalar_map[&ir];
		let se = concatenate(Axis(0), &[sl.slice(s![-2.., ..]), s0.view(), sr.slice(s![..2, ..])])?;

		let gx = ndarray_ops::map_stencil3(&pe, Axis(0), |a, b, c| hydro.plm_gradient(PLM_THETA, a, b, c));
		let gy = ndarray_ops::map_stencil3(&pe, Axis(1), |a, b, c| hydro.plm_gradient(PLM_THETA, a, b, c));
		let gy = ndarray_ops::extend_default_2d(gy, 0, 0, 1, 1);

		let hx = ndarray_ops::map_stencil3(&se, Axis(0), |a, b, c| plm_gradient(PLM_THETA, a, b, c));
		let hy = ndarray_ops::map_stencil3(&se, Axis(1), |a, b, c| plm_gradient(PLM_THETA, a, b, c));
		let hy = ndarray_ops::extend_default_2d(hy, 0, 0, 1, 1);

	    let pxl = pe.slice(s![1..-2, ..]);
	    let pxr = pe.slice(s![2..-1, ..]);
	    let gxl = gx.slice(s![ ..-1, ..]);
	    let gxr = gx.slice(s![1..  , ..]);
	    let sxl = se.slice(s![1..-2, ..]);
	    let sxr = se.slice(s![2..-1, ..]);
	    let hxl = hx.slice(s![ ..-1, ..]);
	    let hxr = hx.slice(s![1..  , ..]);

		let pyl = pe.slice(s![2..-2,  ..-1]);
		let pyr = pe.slice(s![2..-2, 1..  ]);
		let gyl = gy.slice(s![2..-2,  ..-1]);
		let gyr = gy.slice(s![2..-2, 1..  ]);
		let syl = se.slice(s![2..-2,  ..-1]);
		let syr = se.slice(s![2..-2, 1..  ]);
		let hyl = hy.slice(s![2..-2,  ..-1]);
		let hyr = hy.slice(s![2..-2, 1..  ]);

		// let godunov_x = Array::from_shape_fn(pxl.dim(), |i| {
		//     hydro.intercell_flux(pxl[i] + gxl[i] * 0.5, pxr[i] - gxr[i] * 0.5, sxl[i] + hxl[i] * 0.5, sxr[i] - hxr[i] * 0.5, Direction::Radial)
		// });
		// let godunov_y = Array::from_shape_fn(pxl.dim(), |i| {
		//     hydro.intercell_flux(pxl[i] + gxl[i] * 0.5, pxr[i] - gxr[i] * 0.5, sxl[i] + hxl[i] * 0.5, sxr[i] - hxr[i] * 0.5, Direction::Polar)
		// });
	}

	state.time += 0.1;
	state.iteration += 1;

	// std::thread::sleep(std::time::Duration::from_millis(100));

	Ok(())
}
