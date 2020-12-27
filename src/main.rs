/**
 * The Clemson Kilonova Code
 */




// ============================================================================
static DESCRIPTION: &str = env!("CARGO_PKG_DESCRIPTION");
static VERSION_AND_BUILD: &str = git_version::git_version!(prefix=concat!("v", env!("CARGO_PKG_VERSION"), " "));




// ============================================================================
mod mesh;
mod models;
mod physics;
mod scheme;
mod state;
mod tasks;
mod traits;




// ============================================================================
use std::fs::{
    File,
};
use serde::{
    Serialize,
    Deserialize,
};
use mesh::{
    Mesh,
};
use models::{
    JetInCloud,
    HaloKilonova,
};
use physics::{
    AgnosticPrimitive,
    RelativisticHydrodynamics,
};
use state::{
    State,
};
use traits::{
    Conserved,
    Hydrodynamics,
    InitialModel,
};
use tasks::{
    Tasks,
};




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
 * Simulation control: how long to run for, how frequently to perform side
 * effects, etc.
 */
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Control {
    pub final_time: f64,
    pub checkpoint_interval: f64,
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
    pub control: Control,
}


/**
 * App state
 */
#[derive(Serialize, Deserialize)]
struct App {
    state: AgnosticState,
    tasks: Tasks,
    config: Configuration,
}




// ============================================================================
impl App {

    /**
     * Construct a new App instance from a user configuration
     */
    fn from_config(config: Configuration) -> anyhow::Result<Self> {
        let geometry = config.mesh.grid_blocks_geometry()?;
        let state = match &config.hydro {
            AgnosticHydrodynamics::Euler => anyhow::bail!("hydro: euler is not implemented yet"),
            AgnosticHydrodynamics::Relativistic(hydro) => {
                let state = State::from_model(&config.model, hydro, &geometry);
                AgnosticState::Relativistic(state)
            }
        };
        let tasks = Tasks::new();
        Ok(Self{state, tasks, config})
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




// ============================================================================
fn side_effects<C: Conserved>(state: &State<C>, tasks: &mut Tasks, _control: &Control) {

    let mzps = 1e-6 * state.total_zones() as f64 / tasks.iteration_message.lap_seconds();

    println!("[{:05}] t={:.3} blocks={} mzps={:.2})", state.iteration, state.time, state.solution.len(), mzps);

    tasks.write_checkpoint.advance(1.0);
}




// ============================================================================
fn run<H, C>(
    mut state: State<C>,
    mut tasks: Tasks,
    mesh: Mesh,
    hydro: H,
    control: Control) -> anyhow::Result<()>
where
    H: Hydrodynamics<Conserved = C>,
    C: Conserved {
    while state.time < control.final_time {
        scheme::advance(&mut state, &hydro, &mesh);
        side_effects(&state, &mut tasks, &control);
    }
    Ok(())
}




// ============================================================================
fn main() -> anyhow::Result<()> {

    let App{state, tasks, config} = App::build()?;
    let Configuration{hydro, mesh, control, ..} = config;

    println!("{}", DESCRIPTION);
    println!("{}", VERSION_AND_BUILD);

    match (state, hydro) {
        (AgnosticState::Euler, _) => {
            anyhow::bail!("Euler hydrodynamics not implemented")
        },
        (AgnosticState::Relativistic(state), AgnosticHydrodynamics::Relativistic(hydro)) => {
            run(state, tasks, mesh, hydro, control)
        },
        _ => unreachable!(),
    }
}
