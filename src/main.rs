/**
 * The Clemson Kilonova Code
 */


static DESCRIPTION: &str = env!("CARGO_PKG_DESCRIPTION");
static VERSION_AND_BUILD: &str = git_version::git_version!(prefix=concat!("v", env!("CARGO_PKG_VERSION"), " "));


mod mesh;
mod models;
mod physics;
mod products;
mod scheme;
mod state;
mod tasks;
mod traits;


/**
 * Standard and external dependencies
 */
use std::{
    ffi::OsStr,
    fs::{
        File,
        read_to_string,
    },
    io::BufWriter,
    io::BufReader,
    path::Path,
};
use serde::{
    Serialize,
    Deserialize,
};
use enum_dispatch::enum_dispatch;


/**
 * Imports from crate modules
 */
use mesh::Mesh;
use models::{
    JetInCloud,
    HaloKilonova,
};
use physics::{
    AgnosticPrimitive,
    RelativisticHydro,
};
use products::Products;
use state::State;
use traits::{
    Conserved,
    Hydrodynamics,
    InitialModel,
};
use tasks::Tasks;


/**
 * Model choice
 */
#[enum_dispatch(InitialModel)]
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "snake_case")]
pub enum Model {
    JetInCloud(JetInCloud),
    HaloKilonova(HaloKilonova),
}


/**
 * Enum for any of the supported hydrodynamics types
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "snake_case")]
pub enum AgnosticHydro {
    Euler,
    Relativistic(RelativisticHydro),
}


/**
 * Enum for the solution state of any of the supported hydrodynamics types
 */
#[derive(Clone, Serialize, Deserialize)]
pub enum AgnosticState {
    Euler,
    Relativistic(State<hydro_srhd::srhd_2d::Conserved>),
}


/**
 * Simulation control: how long to run for, how frequently to perform side
 * effects, etc
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Control {
    pub final_time: f64,
    pub start_time: f64,
    pub checkpoint_interval: f64,
    pub products_interval: f64,
    pub fold: usize,
    pub num_threads: usize,
    pub snappy_compression: bool,
}


/**
 * User configuration
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Configuration {
    pub hydro: AgnosticHydro,
    pub model: Model,
    pub mesh: Mesh,
    pub control: Control,
}


/**
 * App state
 */
#[derive(Clone, Serialize, Deserialize)]
pub struct App {
    state: AgnosticState,
    tasks: Tasks,
    config: Configuration,
    version: String,
}




// ============================================================================
impl From<State<hydro_srhd::srhd_2d::Conserved>> for AgnosticState {
    fn from(state: State<hydro_srhd::srhd_2d::Conserved>) -> Self {
        Self::Relativistic(state)
    }
}

impl From<RelativisticHydro> for AgnosticHydro {
    fn from(hydro: RelativisticHydro) -> Self {
        Self::Relativistic(hydro)
    }
}

impl AgnosticHydro {
    fn validate(&self) -> anyhow::Result<()> {
        match self {
            AgnosticHydro::Euler => Ok(()),
            AgnosticHydro::Relativistic(hydro) => hydro.validate(),
        }        
    }
}

impl Control {
    fn validate(&self) -> anyhow::Result<()> {
        if self.num_threads == 0 || self.num_threads >= 1024 {
            anyhow::bail!("num_threads must be > 0 and < 1024")
        }
        if self.checkpoint_interval < 0.0 {
            anyhow::bail!("checkpoint_interval <= 0.0")
        }
        if self.products_interval < 0.0 {
            anyhow::bail!("products_interval <= 0.0")
        }
        Ok(())
    }
}




// ============================================================================
impl Configuration {
    fn package<H>(hydro: &H, model: &Model, mesh: &Mesh, control: &Control) -> Self
    where
        H: Hydrodynamics,
        AgnosticHydro: From<H> {
        Configuration{
            hydro: AgnosticHydro::from(hydro.clone()),
            model: model.clone(),
            mesh: mesh.clone(),
            control: control.clone(),
        }
    }

    fn validate(&self) -> anyhow::Result<()> {
        self.hydro.validate()?;
        self.model.validate()?;
        self.mesh.validate(self.control.start_time)?;
        self.control.validate()?;
        Ok(())
    }
}




// ============================================================================
impl App {

    /**
     * Return self as a result, which will be in an error state if any of the
     * configuration items did not pass validation.
     */
    fn validate(self) -> anyhow::Result<Self> {
        self.config.validate()?;
        Ok(self)
    }

    /**
     * Construct a new App instance from a user configuration.
     */
    fn from_config(config: Configuration) -> anyhow::Result<Self> {
        let geometry = config.mesh.grid_blocks_geometry(config.control.start_time);
        let state = match &config.hydro {
            AgnosticHydro::Euler => {
                anyhow::bail!("hydro: euler is not implemented yet")
            },
            AgnosticHydro::Relativistic(hydro) => {
                AgnosticState::from(State::from_model(&config.model, hydro, &geometry, config.control.start_time))
            },
        };
        let tasks = Tasks::new(config.control.start_time);
        Ok(Self{state, tasks, config, version: VERSION_AND_BUILD.to_string()})
    }

    /**
     * Construct a new App instance from a file: may be a config.yaml or a
     * chkpt.0000.cbor.
     */
    fn from_file(filename: &str) -> anyhow::Result<Self> {
        match Path::new(&filename).extension().and_then(OsStr::to_str) {
            Some("yaml") => Self::from_config(serde_yaml::from_str(&read_to_string(filename)?)?),
            Some("cbor") => read_cbor(filename, false),
            Some("cboz") => read_cbor(filename, true),
            _ => anyhow::bail!("unknown input file type '{}'", filename),
        }
    }

    /**
     * Construct a new App instance from a preset (hard-coded) configuration
     * name, or otherwise an input file if no matching preset is found.
     */
    fn from_preset_or_file(input: &str) -> anyhow::Result<Self> {
        match input {
            "jet_in_cloud" => Self::from_config(serde_yaml::from_str(std::include_str!("../setups/jet_in_cloud.yaml"))?),
            _ => Self::from_file(input),
        }
    }

    /**
     * Construct a new App instance from references to the member variables.
     */
    fn package<C, H>(state: &State<C>, tasks: &mut Tasks, hydro: &H, model: &Model, mesh: &Mesh, control: &Control) -> Self
    where
        H: Hydrodynamics<Conserved = C>,
        C: Conserved,
        AgnosticState: From<State<C>>,
        AgnosticHydro: From<H> {
        Self{
            state: AgnosticState::from(state.clone()),
            tasks: tasks.clone(),
            config: Configuration::package(hydro, model, mesh, control),
            version: VERSION_AND_BUILD.to_string(),
        }
    }
}




// ============================================================================
fn parent_directory(path_str: &str) -> String {
    match Path::new(&path_str).parent().and_then(Path::to_str) {
        None     => ".",
        Some("") => ".",
        Some(parent) => parent,
    }.into()
}

fn write_cbor<T: Serialize>(value: &T, path_str: &str, snappy_compression: bool) -> anyhow::Result<()> {
    println!("write {}", path_str);
    let file = File::create(&path_str)?;
    let buffer = BufWriter::new(file);

    if snappy_compression {
        serde_cbor::to_writer(snap::write::FrameEncoder::new(buffer), &value)?;
    } else {
        serde_cbor::to_writer(buffer, &value)?;        
    }
    Ok(())
}

fn read_cbor<T: for<'de> Deserialize<'de>>(path_str: &str, snappy_compression: bool) -> anyhow::Result<T> {
    let file = File::open(path_str)?;
    let buffer = BufReader::new(file);

    if snappy_compression {
        Ok(serde_cbor::from_reader(snap::read::FrameDecoder::new(buffer))?)
    } else {
        Ok(serde_cbor::from_reader(buffer)?)        
    }
}




// ============================================================================
fn side_effects<C, H>(state: &State<C>, tasks: &mut Tasks, hydro: &H, model: &Model, mesh: &Mesh, control: &Control, outdir: &str)
    -> anyhow::Result<()>
where
    H: Hydrodynamics<Conserved = C>,
    C: Conserved,
    AgnosticState: From<State<C>>,
    AgnosticHydro: From<H> {

    let extension = if control.snappy_compression {
        "cboz"
    } else {
        "cbor"
    };

    if tasks.iteration_message.next_time <= state.time {
        let time = tasks.iteration_message.advance(0.0);
        let mzps = 1e-6 * state.total_zones() as f64 / time * control.fold as f64;
        if tasks.iteration_message.count_this_run > 1 {
            println!("[{:05}] t={:.5} blocks={} Mzps={:.2})", state.iteration, state.time, state.solution.len(), mzps);
        }
    }

    if tasks.write_products.next_time <= state.time {
        tasks.write_products.advance(control.products_interval);
        let filename = format!("{}/prods.{:04}.{}", outdir, tasks.write_products.count - 1, extension);
        let config = Configuration::package(hydro, model, mesh, control);
        let products = Products::from_state(state, hydro, mesh, &config)?;
        write_cbor(&products, &filename, control.snappy_compression)?;
    }

    if tasks.write_checkpoint.next_time <= state.time {
        tasks.write_checkpoint.advance(control.checkpoint_interval);
        let filename = format!("{}/chkpt.{:04}.{}", outdir, tasks.write_checkpoint.count - 1, extension);
        let app = App::package(state, tasks, hydro, model, mesh, control);
        write_cbor(&app, &filename, control.snappy_compression)?;
    }

    Ok(())
}




// ============================================================================
fn run<C, H>(mut state: State<C>, mut tasks: Tasks, hydro: H, model: Model, mesh: Mesh, control: Control, outdir: String)
    -> anyhow::Result<()>
where
    H: Hydrodynamics<Conserved = C>,
    C: Conserved,
    AgnosticState: From<State<C>>,
    AgnosticHydro: From<H> {

    let mut block_geometry = mesh.grid_blocks_geometry(state.time);
    let runtime = tokio::runtime::Builder::new_multi_thread()
        .worker_threads(control.num_threads)
        .build()?;

    while state.time < control.final_time {
        side_effects(&state, &mut tasks, &hydro, &model, &mesh, &control, &outdir)?;
        state = scheme::advance(state, &hydro, &model, &mesh, &mut block_geometry, &runtime, control.fold)?;
    }

    side_effects(&state, &mut tasks, &hydro, &model, &mesh, &control, &outdir)?;

    Ok(())
}




// ============================================================================
fn main() -> anyhow::Result<()> {

    let input = match std::env::args().nth(1) {
        None => anyhow::bail!("no input file given"),
        Some(input) => input,
    };
    let outdir = parent_directory(&input);

    println!();
    println!("\t{}", DESCRIPTION);
    println!("\t{}", VERSION_AND_BUILD);
    println!();
    println!("\tinput file ........ {}", input);
    println!("\toutput drectory ... {}", outdir);

    let App{state, tasks, config, ..} = App::from_preset_or_file(&input)?.validate()?;
    let Configuration{hydro, model, mesh, control} = match std::env::args().nth(2) {
        Some(extra) => serde_yaml::from_str(&read_to_string(extra)?)?,
        None => config,
    };

    match (state, hydro) {
        (AgnosticState::Euler, _) => {
            anyhow::bail!("Euler hydrodynamics not implemented")
        },
        (AgnosticState::Relativistic(state), AgnosticHydro::Relativistic(hydro)) => {
            run(state, tasks, hydro, model, mesh, control, outdir)
        },
        _ => unreachable!(),
    }
}
