#![allow(unused)]
use std::fs::File;
use serde::{Serialize, Deserialize};




// ============================================================================
pub static DESCRIPTION: &str = env!("CARGO_PKG_DESCRIPTION");
pub static VERSION_AND_BUILD: &str = git_version::git_version!(prefix=concat!("v", env!("CARGO_PKG_VERSION"), " "));




// ============================================================================
mod mesh;
mod physics;
mod scheme;
mod state;
mod tasks;
mod traits;
mod models;

use std::collections::HashMap;
use mesh::Mesh;
use models::{JetInCloud, HaloKilonova};
use physics::AgnosticPrimitive;
use state::State;
use traits::{Primitive, InitialModel};


/**
 * Model choice
 */
#[enum_dispatch::enum_dispatch(InitialModel)]
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "snake_case")]
enum Model {
    JetInCloud(JetInCloud),
    HaloKilonova(HaloKilonova),
}


/**
 * Mirror of a user configuration
 */
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct UserConfiguration {
    pub model: Model,
    pub mesh: Mesh,
}




// ============================================================================
fn main() -> anyhow::Result<()> {

    let user: UserConfiguration = if let Some(f) = std::env::args().skip(1).next() {
        serde_yaml::from_reader(File::open(f)?)?
    } else {
        anyhow::bail!("no config file given")
    };
    println!("{}\n", serde_yaml::to_string(&user)?.replace("---", ""));

    let system = physics::RelativisticHydrodynamics::new();

    let mut primitive_map = HashMap::new();

    for (index, subgrid) in user.mesh.grid_blocks() {
        println!("{:?} {}", index, subgrid.extent.inner_radius);
        primitive_map.insert(index, physics::grid_primitive(&subgrid, &system, &user.model));
    }

    Ok(())
}
