static DESCRIPTION: &str = env!("CARGO_PKG_DESCRIPTION");
static VERSION_AND_BUILD: &str = git_version::git_version!(prefix=concat!("v", env!("CARGO_PKG_VERSION"), " "));

mod mesh;
mod models;
mod physics;
mod scheme;
mod state;
mod tasks;
mod traits;

use std::fs::File;
use std::collections::HashMap;
use serde::{Serialize, Deserialize};
use mesh::Mesh;
use models::{JetInCloud, HaloKilonova};
use physics::{AgnosticPrimitive, RelativisticHydrodynamics};
use state::State;
use traits::InitialModel;
use tasks::Tasks;




/**
 * Model choice
 */
#[enum_dispatch::enum_dispatch(InitialModel)]
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "snake_case", tag = "type")]
enum Model {
    JetInCloud(JetInCloud),
    HaloKilonova(HaloKilonova),
}


/**
 * Enum for any of the supported hydrodynamics types
 */
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "snake_case", tag = "type")]
enum AgnosticHydrodynamics {
    Euler,
    Relativistic(RelativisticHydrodynamics),
}


/**
 * Enum for the solution state of any of the supported hydrodynamics types
 */
#[derive(Serialize, Deserialize)]
enum AgnosticState {
    Euler,
    Relativistic(State<hydro_srhd::srhd_2d::Conserved>),
}


/**
 * User configuration
 */
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Configuration {
    pub hydro: AgnosticHydrodynamics,
    pub model: Model,
    pub mesh: Mesh,
}


/**
 * App state
 */
#[derive(Serialize, Deserialize)]
struct App {
    state: AgnosticState,
    tasks: Tasks,
    model: Model,
}




// ============================================================================
impl App {

    /**
     * Construct a new App instance from a user configuration
     */
    fn from_config(config: Configuration) -> anyhow::Result<Self> {
        let state = match config.hydro {
            AgnosticHydrodynamics::Euler => anyhow::bail!("hydro: euler is not implemented yet"),
            AgnosticHydrodynamics::Relativistic(_) => AgnosticState::Relativistic(
                State{
                    time: 0.0,
                    iteration: num::rational::Rational64::new(0, 1),
                    solution: HashMap::new(),
                }
            )
        };

        let tasks = Tasks::new();
        let model = config.model;

        Ok(Self{state, tasks, model})
    }

    /**
     * Construct a new App instance from the command line arguments
     */
    fn build() -> anyhow::Result<Self> {
        if let Some(input_file) = std::env::args().skip(1).next() {
            if input_file.ends_with(".yaml") {
                Self::from_config(serde_yaml::from_reader(File::open(input_file)?)?)
            } else if input_file.ends_with(".pk") {
                Ok(serde_pickle::from_reader(File::open(input_file)?)?)
            } else {
                anyhow::bail!("unknown input file type '{}'", input_file)
            }
        }
        else {
            anyhow::bail!("no input file given")
        }
    }
}




/**
 * Perform side effects and mutate the tasks struct
 */
fn side_effects(tasks: &mut Tasks) {
    tasks.write_checkpoint.advance(1.0);    
}




/**
 * Main function
 */
fn main() -> anyhow::Result<()> {

    let App{mut state, mut tasks, model} = App::build()?;

    println!("{}", DESCRIPTION);
    println!("{}", VERSION_AND_BUILD);

    // let mut primitive_map = HashMap::new();

    // for (index, subgrid) in user.mesh.grid_blocks() {
    //     println!("{:?} {}", index, subgrid.extent.inner_radius);
    //     primitive_map.insert(index, physics::grid_primitive(&subgrid, &system, &user.model));
    // }

    for _ in 0..100 {
        side_effects(&mut tasks);
    }

    Ok(())
}
