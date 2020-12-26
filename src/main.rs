#![allow(unused)]
use std::collections::HashMap;
use std::fs::File;
use std::default::Default;
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

use mesh::Mesh;


/**
 * Jet propagating through a kilonova debris cloud and surrounding relativistic
 * envelop
 */
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct JetInCloud {
    pub engine_duration: f64,
    pub engine_strength: f64,
}


/**
 * Explosion in a horizontally stratified external medium
 */
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct HaloKilonova {
    pub explosion_energy: f64,
}


/**
 * Model choice
 */
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
#[serde(rename_all = "snake_case")]
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
    println!("{}", serde_yaml::to_string(&user)?);


    Ok(())
}
