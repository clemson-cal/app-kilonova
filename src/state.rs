use std::collections::HashMap;
use num::ToPrimitive;
use num::rational::Rational64;
use serde::{Serialize, Deserialize};
use ndarray::{Array, ArcArray, Ix2};
use godunov_core::runge_kutta;
use crate::physics::HydroError;
use crate::traits::{
    Conserved,
    Hydrodynamics,
    InitialModel,
    Primitive,
};
use crate::mesh::{
    BlockIndex,
    Mesh,
    GridGeometry
};




/**
 * The solution state for an individual grid block
 */
#[derive(Clone, Serialize, Deserialize)]
pub struct BlockState<C: Conserved> {
    pub conserved: ArcArray<C, Ix2>,
    pub scalar_mass: ArcArray<f64, Ix2>,
}




/**
 * The full solution state for the simulation
 */
#[derive(Clone, Serialize, Deserialize)]
pub struct State<C: Conserved> {
    pub time: f64,
    pub iteration: Rational64,
    pub solution: HashMap<BlockIndex, BlockState<C>>,
}




// ============================================================================
impl<C: Conserved> BlockState<C> {

    /**
     * Generate a block state from the given initial model, hydrodynamics
     * instance and grid geometry.
     */
    pub fn from_model<M, H>(model: &M, hydro: &H, geometry: &GridGeometry, time: f64) -> Self
    where
        M: InitialModel,
        H: Hydrodynamics<Conserved = C>
    {
        let scalar      = geometry.cell_centers.mapv(|c| model.scalar_at(c, time));
        let primitive   = geometry.cell_centers.mapv(|c| hydro.interpret(&model.primitive_at(c, time)));
        let conserved   = primitive.mapv(|p| hydro.to_conserved(p)) * &geometry.cell_volumes;
        let scalar_mass = conserved.mapv(|u| u.lab_frame_mass()) * scalar;

        Self {
            conserved: conserved.to_shared(),
            scalar_mass: scalar_mass.to_shared()
        }
    }

    /**
     * Try to convert the array of conserved quantities in this block to an
     * array of primitive quantities, and return an error if the conversion
     * failed anyhere. This function will not panic.
     */
    pub fn try_to_primitive<H, P>(
        &self,
        hydro: &H,
        geometry: &GridGeometry) -> anyhow::Result<Array<P, Ix2>, HydroError>
    where
        H: Hydrodynamics<Conserved = C, Primitive = P>,
        C: Conserved,
        P: Primitive  
    {
        let u = &self.conserved / &geometry.cell_volumes;
        let x: Result<Vec<_>, _> = u
            .iter()
            .zip(geometry.cell_centers.iter())
            .map(|(&u, &rq)| hydro
                .try_to_primitive(u)
                .map_err(|e| e.at_position(rq)))
            .collect();
        Ok(Array::from_shape_vec(u.dim(), x?).unwrap())
    }
}




// ============================================================================
impl<C: Conserved> State<C> {

    /**
     * Generate a state from the given initial model, hydrodynamics instance,
     * and map of grid geometry.
     */
    pub fn from_model<M, H>(model: &M, hydro: &H, geometry: &HashMap<BlockIndex, GridGeometry>, time: f64) -> Self
    where
        M: InitialModel,
        H: Hydrodynamics<Conserved = C>
    {
        let iteration = Rational64::new(0, 1);
        let solution = geometry.iter().map(|(&i, g)| (i, BlockState::from_model(model, hydro, g, time))).collect();
        Self{time, iteration, solution}
    }

    /**
     * Return the total number of grid zones in this state.
     */
    pub fn total_zones(&self) -> usize {
        self.solution.values().map(|solution| solution.conserved.len()).sum()
    }

    /**
     * Return the indexes of "ghost blocks" just inside and outside the mesh
     * radial extent.
     */
    pub fn inner_outer_boundary_indexes(&self) -> (BlockIndex, BlockIndex) {
        self.min_max_block_indexes_offset_by(1)
    }

    /**
     * Return the indexes of the innermost and outermost block indexes.
     */
    pub fn inner_outer_block_indexes(&self) -> (BlockIndex, BlockIndex) {
        self.min_max_block_indexes_offset_by(0)
    }

    /**
     * Return the time step size, computed from the mesh, the hydrodynamics
     * state, and internal parameters such as the CFL number.
     */
    pub fn time_step<H>(&self, hydro: &H, mesh: &Mesh) -> Result<f64, HydroError>
    where
        H: Hydrodynamics<Conserved = C>
    {
        if let Some(max_signal_speed) = hydro.global_signal_speed() {
            let (index, ..) = self.inner_outer_block_indexes();
            Ok(hydro.cfl_number() * mesh.smallest_spacing(index) / max_signal_speed)
        } else {
            Ok(self.solution.iter().try_fold(f64::MAX, |dt, (index, state)| {
                let geometry = mesh.subgrid(*index).geometry();
                let block_dt = state
                    .try_to_primitive(hydro, &geometry)?
                    .iter()
                    .zip(&geometry.cell_linear_dimension())
                    .fold(dt, |dt, (p, dl)| dt.min(dl / hydro.max_signal_speed(*p))
                );
                Ok(dt.min(block_dt))
            })? * hydro.cfl_number())
        }
    }

    fn min_max_block_indexes_offset_by(&self, delta: i32) -> (BlockIndex, BlockIndex) {
        let mut min = (i32::MAX, 0);
        let mut max = (i32::MIN, 0);
        for i in self.solution.keys() {
            min = (min.0.min(i.0 - delta), min.1);
            max = (max.0.max(i.0 + delta), max.1);
        }
        (min, max)
    }
}




// ============================================================================
impl<C: Conserved> runge_kutta::WeightedAverage for BlockState<C> {
    fn weighted_average(self, br: Rational64, s0: &Self) -> Self {
        let s1 = self;
        let bf = br.to_f64().unwrap();
        let u0 = s0.conserved.clone();
        let u1 = s1.conserved.clone();
        let c0 = s0.scalar_mass.clone();
        let c1 = s1.scalar_mass.clone();

        Self {
            conserved:   u1 * (-bf + 1.) + u0 * bf,
            scalar_mass: c1 * (-bf + 1.) + c0 * bf,
        }
    }
}




// ============================================================================
impl<C: Conserved> runge_kutta::WeightedAverage for State<C> {
    fn weighted_average(self, br: Rational64, s0: &Self) -> Self {
        let bf = br.to_f64().unwrap();
        let s_avg = self.solution
            .into_iter()
            .map(|(index, s1)| (index, s1.weighted_average(br, &s0.solution[&index])));

        Self {
            time:      self.time      * (-bf + 1.) + s0.time      * bf,
            iteration: self.iteration * (-br + 1 ) + s0.iteration * br,
            solution: s_avg.into_iter().collect(),
        }
    }
}




// ============================================================================
#[async_trait::async_trait]
impl<C: Conserved> runge_kutta::WeightedAverageAsync for State<C> {

    type Runtime = tokio::runtime::Runtime;

    async fn weighted_average(self, br: Rational64, s0: &Self, runtime: &Self::Runtime) -> Self {
        use futures::future::join_all;
        use godunov_core::runge_kutta::WeightedAverage;

        let bf = br.to_f64().unwrap();
        let s_avg = self.solution.into_iter().map(|(index, s1)| {
            let s0 = s0.clone();
            async move {
                runtime.spawn(
                    async move {
                        (index, s1.weighted_average(br, &s0.solution[&index]))
                    }
                ).await.unwrap()
            }
        });

        Self {
            time:      self.time      * (-bf + 1.) + s0.time      * bf,
            iteration: self.iteration * (-br + 1 ) + s0.iteration * br,
            solution: join_all(s_avg).await.into_iter().collect(),
        }
    }
}
