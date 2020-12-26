use std::collections::HashMap;
use ndarray::{ArcArray, Ix1, Ix2};
use serde::{Serialize, Deserialize};


/**
 * Type alias for a 2D block index
 */
pub type BlockIndex = (usize, usize);


/**
 * A struct to hold cached calculations of geometric primitives relevant to
 * spherical polar geometry
 */
#[derive(Clone)]
pub struct BlockGeometry {
    radial_vertices:   ArcArray<f64, Ix1>,
    polar_vertices:    ArcArray<f64, Ix1>,
    cell_centers:      ArcArray<(f64, f64), Ix2>,
    cell_volumes:      ArcArray<(f64, f64), Ix2>,
    radial_face_areas: ArcArray<(f64, f64), Ix2>,
    polar_face_areas:  ArcArray<(f64, f64), Ix2>,
}


/**
 * A volume in (r, theta) space; polar section of a spherical annulus
 */
#[derive(Clone)]
pub struct SphericalPolarExtent {
    pub inner_radius: f64,
    pub outer_radius: f64,
    pub lower_theta: f64,
    pub upper_theta: f64,
}


/**
 * A spherical polar extent, together with zone counts to make even subdivisions
 * in log-r and polar angle
 */
#[derive(Clone)]
pub struct SphericalPolarGrid {
    pub extent: SphericalPolarExtent,
    pub num_zones_r: usize,
    pub num_zones_q: usize,
}


/**
 * Abstract description of a spherical polar mesh
 */
#[derive(Serialize, Deserialize)] #[serde(deny_unknown_fields)]
pub struct Mesh {

    /// Inner boundary radius; the grid will start precisely here
    pub inner_radius: f64,

    /// Outer radius; the grid may extend up to one block past this
    pub outer_radius: f64,

    /// Number of zones from pole to pole
    pub num_polar_zones: usize,

    /// Number of radial zones in each block
    pub block_size: usize,
}




// ============================================================================
impl SphericalPolarExtent {
    pub fn grid(&self, num_zones_r: usize, num_zones_q: usize) -> SphericalPolarGrid {
        SphericalPolarGrid{
            extent: self.clone(),
            num_zones_r,
            num_zones_q,
        }
    }
}




// ============================================================================
impl SphericalPolarGrid {
    fn vertex_coordinate(&self, i: i64, j: i64) -> (f64, f64) {
        let (y0, y1) = (self.extent.inner_radius.log(10.0), self.extent.outer_radius.log(10.0));
        let (q0, q1) = (self.extent.lower_theta, self.extent.upper_theta);
        let dy = (y1 - y0) / self.num_zones_r as f64;
        let dq = (q1 - q0) / self.num_zones_q as f64;
        let y = y0 + dy * i as f64;
        let q = q0 + dq * j as f64;
        (y.powf(10.0), q)
    }

    fn zone(&self, i: i64, j: i64) -> SphericalPolarExtent {
        let lower = self.vertex_coordinate(i + 0, j + 0);
        let upper = self.vertex_coordinate(i + 1, j + 1);
        SphericalPolarExtent{
            inner_radius: lower.0,
            outer_radius: upper.0,
            lower_theta: lower.1,
            upper_theta: upper.1,
        }
    }
}




// ============================================================================
impl Mesh {
    pub fn grid_blocks(&self) -> HashMap<BlockIndex, SphericalPolarGrid> {
        let mut i = 0;
        let mut r = self.inner_radius;
        let mut blocks = HashMap::new();
        let block_dlogr = 10.0;

        todo!("calculate block_dlogr");

        while r < self.outer_radius {
            let extent = SphericalPolarExtent{
                inner_radius: r,
                outer_radius: r * block_dlogr,
                lower_theta: 0.0,
                upper_theta: std::f64::consts::PI,
            };
            blocks.insert((i, 0), extent.grid(self.block_size, self.num_polar_zones));

            r *= block_dlogr;
            i += 1;
        }
        return blocks;
    }
}
