use std::collections::HashMap;
use ndarray::{Array, Axis, concatenate, s};
use tokio::runtime::Runtime;
use crate::mesh::{BlockIndex, GridGeometry, Mesh};
use crate::Model;
use crate::physics::Direction;
use crate::state::{State, BlockState};
use crate::traits::{Conserved, Hydrodynamics};




// ============================================================================
pub fn advance_rk<H, C>(state: State<C>, hydro: &H, model: &Model, mesh: &Mesh, geometry: &HashMap<BlockIndex, GridGeometry>, dt: f64, runtime: &Runtime)
    -> anyhow::Result<State<C>>
where
    H: Hydrodynamics<Conserved = C>,
    C: Conserved {

    use futures::FutureExt;
    use futures::future::join_all;

    let mut stage_map = HashMap::new();
    let mut new_state_vec = Vec::new();
    let mut stage_primitive_and_scalar = |index: BlockIndex, state: BlockState<C>, hydro: H, geometry: GridGeometry| {
        let stage = async move {
            let s = (state.scalar_mass / state.conserved.map(Conserved::lab_frame_mass)).to_shared();
            let p = (state.conserved / geometry.cell_volumes).mapv(|q| hydro.to_primitive(q)).to_shared();
            (p, s)
        };
        stage_map.insert(index, runtime.spawn(stage).map(|f| f.unwrap()).shared());
    };


    for (index, state) in &state.solution {
        stage_primitive_and_scalar(index.clone(), state.clone(), hydro.clone(), geometry[index].clone())
    }


    let (inner_bnd_index, outer_bnd_index) = state.inner_outer_boundary_indexes();
    let inner_bnd_geom = mesh.subgrid(inner_bnd_index).geometry();
    let outer_bnd_geom = mesh.subgrid(outer_bnd_index).geometry();
    let inner_bnd_state = BlockState::from_model(model, hydro, &inner_bnd_geom, state.time);
    let outer_bnd_state = BlockState::from_model(model, hydro, &outer_bnd_geom, state.time);
    stage_primitive_and_scalar(inner_bnd_index, inner_bnd_state, hydro.clone(), inner_bnd_geom);
    stage_primitive_and_scalar(outer_bnd_index, outer_bnd_state, hydro.clone(), outer_bnd_geom);


    for (&index, state) in &state.solution {

        let hydro = hydro.clone();
        let state = state.clone();
        let stage_map = stage_map.clone();
        let geometry = geometry[&index].clone();

        let entry = async move {
            let il = (index.0 - 1, index.1);
            let i0 = (index.0,     index.1);
            let ir = (index.0 + 1, index.1);

            let (pl, sl) = stage_map[&il].clone().await;
            let (p0, s0) = stage_map[&i0].clone().await;
            let (pr, sr) = stage_map[&ir].clone().await;
            let pe = concatenate(Axis(0), &[pl.slice(s![-2.., ..]), p0.view(), pr.slice(s![..2, ..])]).unwrap();
            let se = concatenate(Axis(0), &[sl.slice(s![-2.., ..]), s0.view(), sr.slice(s![..2, ..])]).unwrap();

            let gx = ndarray_ops::map_stencil3(&pe, Axis(0), |a, b, c| hydro.plm_gradient_primitive(a, b, c));
            let gy = ndarray_ops::map_stencil3(&pe, Axis(1), |a, b, c| hydro.plm_gradient_primitive(a, b, c));
            let gy = ndarray_ops::extend_default_2d(gy, 0, 0, 1, 1);

            let hx = ndarray_ops::map_stencil3(&se, Axis(0), |a, b, c| hydro.plm_gradient_scalar(a, b, c));
            let hy = ndarray_ops::map_stencil3(&se, Axis(1), |a, b, c| hydro.plm_gradient_scalar(a, b, c));
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

            let godunov_x = Array::from_shape_fn(pxl.dim(), |i| {
                hydro.intercell_flux(
                    pxl[i] + gxl[i] * 0.5, pxr[i] - gxr[i] * 0.5,
                    sxl[i] + hxl[i] * 0.5, sxr[i] - hxr[i] * 0.5, Direction::Radial)
            });
            let godunov_y = Array::from_shape_fn(pyl.dim(), |i| {
                hydro.intercell_flux(
                    pyl[i] + gyl[i] * 0.5, pyr[i] - gyr[i] * 0.5,
                    syl[i] + hyl[i] * 0.5, syr[i] - hyr[i] * 0.5, Direction::Polar)
            });

            let fx = godunov_x.mapv(|(f, _)| f) * &geometry.radial_face_areas;
            let gx = godunov_x.mapv(|(_, g)| g) * &geometry.radial_face_areas;
            let fy = ndarray_ops::extend_default_2d(godunov_y.mapv(|(f, _)| f), 0, 0, 1, 1) * &geometry.polar_face_areas;
            let gy = ndarray_ops::extend_default_2d(godunov_y.mapv(|(_, g)| g), 0, 0, 1, 1) * &geometry.polar_face_areas;

            let sc = ndarray::azip![
                &p0,
                &geometry.cell_centers,
                &geometry.cell_volumes]
            .apply_collect(|&p, &c, &dv| hydro.geometrical_source_terms(p, c) * dv);

            let du = ndarray::azip![
                &sc,
                fx.slice(s![..-1,..]),
                fx.slice(s![ 1..,..]),
                fy.slice(s![..,..-1]),
                fy.slice(s![.., 1..])]
            .apply_collect(|&s, &a, &b, &c, &d| (s - (b - a) - (d - c)) * dt);

            let ds = ndarray::azip![
                gx.slice(s![..-1,..]),
                gx.slice(s![ 1..,..]),
                gy.slice(s![..,..-1]),
                gy.slice(s![.., 1..])]
            .apply_collect(|&a, &b, &c, &d| ((b - a) + (d - c)) * -dt);

            let new_state = BlockState{
                conserved: (&state.conserved + &du).to_shared(),
                scalar_mass: (&state.scalar_mass + &ds).to_shared(),
            };
            (index, new_state)
        };
        new_state_vec.push(runtime.spawn(entry).map(|f| f.unwrap()));
    }


    Ok(State{
        time: state.time + dt,
        iteration: state.iteration + 1,
        solution: runtime.block_on(join_all(new_state_vec)).into_iter().collect(),
    })
}




// ============================================================================
fn add_remove_blocks<H, C>(state: &mut State<C>, hydro: &H, model: &Model, mesh: &Mesh, geometry: &mut HashMap<BlockIndex, GridGeometry>)
where
    H: Hydrodynamics<Conserved = C>,
    C: Conserved {

    let (inner_index, outer_index) = state.inner_outer_block_indexes();
    let solution = &mut state.solution;

    if mesh.subgrid_extent(inner_index).outer_radius < mesh.inner_excision_surface(state.time) {
        geometry.remove(&inner_index);
        solution.remove(&inner_index);
    }

    if mesh.subgrid_extent(outer_index).outer_radius < mesh.outer_excision_surface(state.time) {
        let new_block_index = (outer_index.0 + 1, outer_index.1);
        let new_block_geometry = mesh.subgrid(new_block_index).geometry();
        let new_block_state = BlockState::from_model(model, hydro, &new_block_geometry, state.time);

        geometry.insert(new_block_index, new_block_geometry);
        solution.insert(new_block_index, new_block_state);
    }    
}




// ============================================================================
pub fn advance<H, C>(mut state: State<C>, hydro: &H, model: &Model, mesh: &Mesh, geometry: &mut HashMap<BlockIndex, GridGeometry>, runtime: &Runtime)
    -> anyhow::Result<State<C>>
where
    H: Hydrodynamics<Conserved = C>,
    C: Conserved {

    let dt = hydro.time_step(&state, mesh);

    if mesh.moving_excision_surfaces() {
        add_remove_blocks(&mut state, hydro, model, mesh, geometry);
    }

    let update = |state| advance_rk(state, hydro, model, mesh, geometry, dt, &runtime).unwrap();
    let runge_kutta = hydro.runge_kutta_order();
    let fold = 1;

    for _ in 0..fold {
        state = runge_kutta.advance(state, update);
    }
    Ok(state)
}
