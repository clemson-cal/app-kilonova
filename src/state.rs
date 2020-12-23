use std::collections::HashMap;
use num::ToPrimitive;
use num::rational::Rational64;
use ndarray::{ArcArray, Ix2};
use godunov_core::runge_kutta;
use crate::traits::Conserved;




// ============================================================================
#[derive(Copy, Clone, Eq, Hash)]
struct BlockIndex {
    i: usize,
    j: usize,
}




// ============================================================================
impl std::cmp::PartialEq for BlockIndex {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i && self.j == other.j
    }
}




// ============================================================================
#[derive(Clone)]
struct BlockState<C: Conserved> {
    conserved: ArcArray<C, Ix2>,
    scalar: ArcArray<f64, Ix2>,
}




// ============================================================================
#[derive(Clone)]
struct State<C: Conserved> {
    time: f64,
    iteration: Rational64,
    solution: HashMap<BlockIndex, BlockState<C>>,
}




// ============================================================================
impl<C: Conserved> runge_kutta::WeightedAverage for BlockState<C>
{
    fn weighted_average(self, br: Rational64, s0: &Self) -> Self
    {
        let s1 = self;
        let bf = br.to_f64().unwrap();
        let u0 = s0.conserved.clone();
        let u1 = s1.conserved.clone();
        let c0 = s0.scalar.clone();
        let c1 = s1.scalar.clone();

        Self{
            conserved: u1 * (-bf + 1.) + u0 * bf,
            scalar:    c1 * (-bf + 1.) + c0 * bf,
        }
    }
}




// ============================================================================
#[async_trait::async_trait]
impl<C: Conserved> runge_kutta::WeightedAverageAsync for State<C>
{
    type Runtime = tokio::runtime::Runtime;
    async fn weighted_average(self, br: Rational64, s0: &Self, runtime: &Self::Runtime) -> Self
    {
        use futures::future::FutureExt;
        use futures::future::join_all;
        use godunov_core::runge_kutta::WeightedAverage;

        let bf = br.to_f64().unwrap();
        let s_avg = self.solution.into_iter().map(|(index, s1)| {
            let s0 = s0.clone();
            runtime.spawn(
                async move {
                    (index, s1.weighted_average(br, &s0.solution[&index]))
                }
            ).map(|f| f.unwrap())
        });

        State{
            time:      self.time      * (-bf + 1.) + s0.time      * bf,
            iteration: self.iteration * (-br + 1 ) + s0.iteration * br,
            solution: join_all(s_avg).await.into_iter().collect(),
        }
    }
}
