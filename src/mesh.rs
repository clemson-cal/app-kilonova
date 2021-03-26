use std::collections::HashMap;
use std::f64::consts::PI;
use ndarray::{ArcArray, Array, Ix1, Ix2};
use serde::{Serialize, Deserialize};




/**
 * Type alias for a 2D block index
 */
pub type BlockIndex = (i32, usize);




/**
 * A struct to hold cached calculations of geometric primitives relevant to
 * spherical polar geometry
 */
#[derive(Clone)]
pub struct GridGeometry {
    pub radial_vertices:   ArcArray<f64, Ix1>,
    pub radial_face_areas: ArcArray<f64, Ix2>,
    pub polar_vertices:    ArcArray<f64, Ix1>,
    pub polar_face_areas:  ArcArray<f64, Ix2>,
    pub cell_volumes:      ArcArray<f64, Ix2>,
    pub cell_centers:      ArcArray<(f64, f64), Ix2>,
}




// ============================================================================
impl GridGeometry {

    /**
     * Return a 2D array of the smallest linear dimension of each grid cell.
     */
    pub fn cell_linear_dimension(&self) -> Array<f64, Ix2> {
        Array::from_shape_fn(self.cell_centers.dim(), |(i, j)| {
            let dr = self.radial_vertices[i + 1] - self.radial_vertices[i];
            let dq = self.polar_vertices[j + 1] - self.polar_vertices[j];
            dr.min(dq * self.radial_vertices[i])
        })
    }
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
#[derive(Clone, Serialize, Deserialize)]
pub struct Mesh {

    /// Radius of the inner surface of the i=0 grid block
    pub reference_radius: f64,

    /// Inner boundary radius; the grid will start precisely here
    pub inner_radius: f64,

    /// Outer radius; the grid may extend up to one block past this
    pub outer_radius: f64,

    /// Speed of the inner excision surface (IES)
    pub inner_excision_speed: f64,

    /// Speed of the outer excision surface (OES)
    pub outer_excision_speed: f64,

    /// Number of radial zones per decade. Use None to select a value that
    /// makes the zones square.
    pub num_radial_zones: Option<usize>,

    /// Number of zones from pole to pole
    pub num_polar_zones: usize,

    /// Number of radial zones in each block
    pub block_size: usize,

    /// Time after which the mesh excision starts
    pub excision_delay: Option<f64>,
}




// ============================================================================
fn cell_volume(c0: (f64, f64), c1: (f64, f64)) -> f64
{
    let dcost = -(f64::cos(c1.1) - f64::cos(c0.1));
    2.0 * PI * (c1.0.powi(3) - c0.0.powi(3)) / 3.0 * dcost
}

fn face_area(c0: (f64, f64), c1: (f64, f64)) -> f64
{
    let s0 = c0.0 * f64::sin(c0.1);
    let s1 = c1.0 * f64::sin(c1.1);
    let z0 = c0.0 * f64::cos(c0.1);
    let z1 = c1.0 * f64::cos(c1.1);
    let ds = s1 - s0;
    let dz = z1 - z0;
    PI * (s0 + s1) * (ds * ds + dz * dz).sqrt()
}




// ============================================================================
impl SphericalPolarExtent {

    /**
     * Create a grid from this r-theta area with the given number of zones in
     * the radial and polar directions.
     */
    pub fn grid(&self, num_zones_r: usize, num_zones_q: usize) -> SphericalPolarGrid {
        SphericalPolarGrid{
            extent: self.clone(),
            num_zones_r,
            num_zones_q,
        }
    }

    /**
     * Return the geometric centroid of this r-theta area.
     */
    pub fn centroid(&self) -> (f64, f64) {
        let q = 0.5 * (self.lower_theta + self.upper_theta);
        let r = (self.inner_radius * self.outer_radius).sqrt();
        (r, q)
    }

    /**
     * Return the 3D volume of this extent.
     */
    pub fn volume(&self) -> f64 {
        let c0 = (self.inner_radius, self.lower_theta);
        let c1 = (self.outer_radius, self.upper_theta);
        cell_volume(c0, c1)
    }

    /**
     * Return the area of the inner radial surface.
     */
    pub fn face_area_r(&self) -> f64 {
        let c0 = (self.inner_radius, self.lower_theta);
        let c1 = (self.inner_radius, self.upper_theta);
        face_area(c0, c1)
    }

    /**
     * Return the area of the lower theta surface.
     */
    pub fn face_area_q(&self) -> f64 {
        let c0 = (self.inner_radius, self.lower_theta);
        let c1 = (self.outer_radius, self.lower_theta);
        face_area(c0, c1)
    }
}




// ============================================================================
impl SphericalPolarGrid {

    /**
     * Return the r-theta vertex coordinate for index (i, j). For this function,
     * (i, j) are allowed to be outside the formal extent of this grid,
     * including being negative.
     */
    pub fn vertex_coordinate_signed(&self, i: i64, j: i64) -> (f64, f64) {
        let (y0, y1) = (self.extent.inner_radius.log(10.0), self.extent.outer_radius.log(10.0));
        let (q0, q1) = (self.extent.lower_theta, self.extent.upper_theta);
        let dy = (y1 - y0) / self.num_zones_r as f64;
        let dq = (q1 - q0) / self.num_zones_q as f64;
        let y = y0 + dy * i as f64;
        let q = q0 + dq * j as f64;
        (f64::powf(10.0, y), q)
    }

    /**
     * Return the r-theta coordinates of the vertex at the given index (i, j).
     */
    pub fn vertex_coordinate(&self, i: usize, j: usize) -> (f64, f64) {
        self.vertex_coordinate_signed(i as i64, j as i64)
    }

    /**
     * Return the spherical polar extent of a cell with index (i, j). For this
     * function, (i, j) are allowed to be outside the formal extent of this
     * grid, including being negative.
     */
    pub fn zone_signed(&self, i: i64, j: i64) -> SphericalPolarExtent {
        let lower = self.vertex_coordinate_signed(i + 0, j + 0);
        let upper = self.vertex_coordinate_signed(i + 1, j + 1);
        SphericalPolarExtent{
            inner_radius: lower.0,
            outer_radius: upper.0,
            lower_theta: lower.1,
            upper_theta: upper.1,
        }
    }

    /**
     * Return the spherical polar extent of a cell with the given index (i, j).
     */
    pub fn zone(&self, index: (usize, usize)) -> SphericalPolarExtent {
        self.zone_signed(index.0 as i64, index.1 as i64)
    }

    /**
     * Return a grid [`GridGeometry`] instance for this grid patch, containing
     * cached geometry data.
     */
    pub fn geometry(&self) -> GridGeometry {
        let nr = self.num_zones_r;
        let nq = self.num_zones_q;
        let radial_vertices   = ArcArray::from_shape_fn(nr + 1, |i| self.vertex_coordinate(i, 0).0);
        let polar_vertices    = ArcArray::from_shape_fn(nq + 1, |j| self.vertex_coordinate(0, j).1);
        let radial_face_areas = ArcArray::from_shape_fn((nr + 1, nq), |index| self.zone(index).face_area_r());
        let polar_face_areas  = ArcArray::from_shape_fn((nr, nq + 1), |index| self.zone(index).face_area_q());
        let cell_volumes      = ArcArray::from_shape_fn((nr, nq), |index| self.zone(index).volume());
        let cell_centers      = ArcArray::from_shape_fn((nr, nq), |index| self.zone(index).centroid());

        GridGeometry{
            radial_vertices,
            radial_face_areas,
            polar_vertices,
            polar_face_areas,
            cell_volumes,
            cell_centers,
        }
    }
}




// ============================================================================
impl Mesh {

    pub fn validate(&self, time: f64) -> anyhow::Result<()> {
        if self.reference_radius <= 0.0 || self.inner_radius < 0.0 || self.outer_radius < 0.0 {
            anyhow::bail!("all radii must be positive")
        }
        if self.reference_radius > self.outer_excision_surface(time) {
            anyhow::bail!("the reference radius is outside the outer excision surface")
        }
        if self.excision_delay.unwrap_or(0.0) < 0.0 {
            anyhow::bail!("the excision delay time must be non-negative")
        }
        if self.inner_excision_speed < 0.0 || self.outer_excision_speed < 0.0 {
            anyhow::bail!("the excision surface speeds must be non-negative")
        }
        if self.outer_excision_speed < self.inner_excision_speed {
            anyhow::bail!("outer_excision_speed < inner_excision_speed (the IES would eventually overtake the OES)")
        }
        if self.block_size < 2 {
            anyhow::bail!("must have at least 2 radial zones per block")
        }
        if self.num_polar_zones != 1 && self.num_polar_zones < 16 {
            anyhow::bail!("must have at least 2 radial zones per block")
        }
        if self.num_polar_zones == 1 && self.num_radial_zones.is_none() {
            anyhow::bail!("num_radial_zones is not optional when num_polar_zones=1")            
        }
        Ok(())
    }

    /**
     * Return the smallest grid spacing on the given block.
     */
    pub fn smallest_spacing(&self, index: BlockIndex) -> f64 {
        let zone = self.subgrid(index).zone((0, 0));
        zone.outer_radius - zone.inner_radius
    }

    /**
     * Return true if either of the IES or the OES have non-zero speeds.
     */
    pub fn moving_excision_surfaces(&self) -> bool {
        self.inner_excision_speed > 0.0 || self.outer_excision_speed > 0.0
    }

    /**
     * Radius of the inner excision surface (IES). The IES is at the
     * `inner_radius` at t=0, and moves outwards at the speed
     * `inner_excision_speed`. Mesh blocks are removed from the mesh if they are
     * fully within the IES.
     */
    pub fn inner_excision_surface(&self, time: f64) -> f64 {
        let t_start = self.excision_delay.unwrap_or(0.0);
        self.inner_radius + (time - t_start).max(0.0) * self.inner_excision_speed
    }

    /**
     * Radius of the outer excision surface (OES). The OES is at the
     * `outer_radius` at t=0, and moves outwards at the speed
     * `outer_excision_speed`. Mesh blocks are added to the mesh if they are
     * fully within by the OES, but not fully within the IES.
     */
    pub fn outer_excision_surface(&self, time: f64) -> f64 {
        let t_start = self.excision_delay.unwrap_or(0.0);
        self.outer_radius + (time - t_start).max(0.0) * self.outer_excision_speed
    }

    /**
     * Return the radial zone spacing, dlogr = log(r1 / r0).
     */
    pub fn zone_dlogr(&self) -> f64 {
        match self.num_radial_zones {
            Some(nr) => 1.0 / nr as f64,
            None => PI / self.num_polar_zones as f64,
        }
    }

    /**
     * Return the number of decades per block.
     */
    pub fn block_dlogr(&self) -> f64 {
        self.block_size as f64 * self.zone_dlogr()
    }

    /**
     * Return the extent of the subgrid at this index.
     */
    pub fn subgrid_extent(&self, index: BlockIndex) -> SphericalPolarExtent {

        let (q0, q1) = if self.num_polar_zones == 1 {
            (PI * 0.5 - self.zone_dlogr(), PI * 0.5 + self.zone_dlogr())
        } else {
            (0.0, PI)
        };

        SphericalPolarExtent {
            inner_radius: self.reference_radius * (1.0 + self.block_dlogr()).powf(index.0 as f64),
            outer_radius: self.reference_radius * (1.0 + self.block_dlogr()).powf(index.0 as f64 + 1.0),
            lower_theta: q0,
            upper_theta: q1,
        }
    }

    /**
     * Return the subgrid object at the given index.
     */
    pub fn subgrid(&self, index: BlockIndex) -> SphericalPolarGrid {
        self.subgrid_extent(index).grid(self.block_size, self.num_polar_zones)
    }

    /**
     * Return a map of the subgrid objects on this mesh.
     */
    pub fn grid_blocks(&self, time: f64) -> HashMap<BlockIndex, SphericalPolarGrid> {
        let mut blocks = HashMap::new();
        for i in 0.. {
            let index = (i, 0);
            let extent = self.subgrid_extent(index);

            if extent.inner_radius >= self.outer_excision_surface(time) {
                break
            } else {
                blocks.insert(index, extent.grid(self.block_size, self.num_polar_zones));
            }
        }
        blocks
    }

    /**
     * Return a map of the subgrid geometry objects on this mesh (for convenience).
     */
    pub fn grid_blocks_geometry(&self, time: f64) -> HashMap<BlockIndex, GridGeometry> {
        self.grid_blocks(time)
            .iter()
            .map(|(&index, grid)| (index, grid.geometry()))
            .collect()
    }
}
