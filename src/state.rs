use std::collections::HashMap;
use num::ToPrimitive;
use num::rational::Rational64;
use serde::{Serialize, Deserialize};
use ndarray::{ArcArray, Ix2};
use godunov_core::runge_kutta;
use crate::traits::{Hydrodynamics, InitialModel, Conserved};
use crate::mesh::{BlockIndex, GridGeometry};




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
    pub fn from_model<M, H>(model: &M, hydro: &H, geometry: &GridGeometry) -> Self
    where
        M: InitialModel,
        H: Hydrodynamics<Conserved = C> {

        let scalar      = geometry.cell_centers.mapv(|c| model.scalar_at(c));
        let primitive   = geometry.cell_centers.mapv(|c| hydro.interpret(&model.primitive_at(c)));
        let conserved   = primitive.mapv(|p| hydro.to_conserved(p)) * &geometry.cell_volumes;
        let scalar_mass = conserved.mapv(|u| u.lab_frame_mass()) * scalar;

        Self{
            conserved: conserved.to_shared(),
            scalar_mass: scalar_mass.to_shared()
        }
    }
}




// ============================================================================
impl<C: Conserved> State<C> {

    /**
     * Generate a state from the given initial model, hydrodynamics instance,
     * and map of grid geometry.
     */
    pub fn from_model<M, H>(model: &M, hydro: &H, geometry: &HashMap<BlockIndex, GridGeometry>) -> Self
    where
        M: InitialModel,
        H: Hydrodynamics<Conserved = C> {

        let time = 0.0;
        let iteration = Rational64::new(0, 1);
        let solution = geometry.iter().map(|(&i, g)| (i, BlockState::from_model(model, hydro, g))).collect();

        Self{time, iteration, solution}
    }

    pub fn total_zones(&self) -> usize {
        self.solution.values().map(|solution| solution.conserved.len()).sum()
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

        Self{
            conserved:   u1 * (-bf + 1.) + u0 * bf,
            scalar_mass: c1 * (-bf + 1.) + c0 * bf,
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

        State{
            time:      self.time      * (-bf + 1.) + s0.time      * bf,
            iteration: self.iteration * (-br + 1 ) + s0.iteration * br,
            solution: join_all(s_avg).await.into_iter().collect(),
        }
    }
}
