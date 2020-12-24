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




// ============================================================================
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Model {
    num_blocks: usize,
    block_size: usize,
}




// ============================================================================
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Execution {
    num_threads: usize,
    outdir: String,
}




// ============================================================================
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Run {
    model: Model,
    execution: Execution,
}




// ============================================================================
fn main() -> anyhow::Result<()> {


    let run: Run = if let Some(config_file) = std::env::args().skip(1).next() {
        serde_yaml::from_reader(File::open(config_file)?)?
    } else {
        anyhow::bail!("no config file given")
    };


    println!("{}", toml::to_string(&run)?);

    // let mut file = std::fs::File::create("chkpt.0000.pk")?;
    // serde_pickle::to_writer(&mut file, &model, true)?;

    Ok(())
}
