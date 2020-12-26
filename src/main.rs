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

    let run: Model = if let Some(config_file) = std::env::args().skip(1).next() {
        serde_yaml::from_reader(File::open(config_file)?)?
    } else {
        anyhow::bail!("no config file given")
    };

    let mut solution = HashMap::new();

    let block_state = state::BlockState{
        conserved: ndarray::ArcArray::default((128, 128)),
        scalar:    ndarray::ArcArray::default((128, 128)),
    };

    for i in 0..96 {
        solution.insert((i, 0), block_state.clone());
    }

    let state = state::State::<hydro_srhd::srhd_2d::Conserved>{
        time: 0.0,
        iteration: num::rational::Rational64::new(0, 1),
        solution: solution,
    };

    println!("{}", toml::to_string(&run)?);

    let file = std::fs::File::create("chkpt.0000.pk")?;
    let mut writer = std::io::BufWriter::new(file);
    println!("serialize with pickle...");
    serde_pickle::to_writer(&mut writer, &state, true)?;
    println!("done!");

    Ok(())
}
