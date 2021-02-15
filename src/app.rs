pub static DESCRIPTION: &str = env!("CARGO_PKG_DESCRIPTION");
pub static VERSION_AND_BUILD: &str = git_version::git_version!(prefix=concat!("v", env!("CARGO_PKG_VERSION"), " "));


use std::{
    ffi::OsStr,
    fs::{File, read_to_string},
    path::Path,
};
use serde::{
    Serialize,
    Deserialize,
};
use yaml_patch::Patch;


use crate::mesh::Mesh;
use crate::models::{
    JetInCloud,
    HaloKilonova,
};
use crate::physics::{
    AnyPrimitive,
    RelativisticHydro,
    NewtonianHydro,
};
use crate::state::State;
use crate::traits::{
    Conserved,
    Hydrodynamics,
    InitialModel,
};
use crate::tasks::Tasks;
use crate::io;


// ============================================================================
#[derive(thiserror::Error, Debug)]
pub enum Error {

    #[error("{0}")]
    IO(#[from] std::io::Error),

    #[error("{0}")]
    SerdeYaml(#[from] serde_yaml::Error),

    #[error("{0}")]
    YamlPatch(#[from] yaml_patch::Error),

    #[error("{0}")]
    AppIO(#[from] io::Error),

    #[error("unknown input file type '{0}'")]
    UnknownInputType(String),
}


/**
 * Model choice
 */
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
pub enum AnyHydro {
    Newtonian(NewtonianHydro),
    Relativistic(RelativisticHydro),
}


/**
 * Enum for the solution state of any of the supported hydrodynamics types
 */
#[derive(Clone, Serialize, Deserialize)]
pub enum AnyState {
    Newtonian(State<hydro_euler::euler_2d::Conserved>),
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

    /// Deprecated
    #[serde(default)]
    pub snappy_compression: bool,
}


/**
 * User configuration
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Configuration {
    pub hydro: AnyHydro,
    pub model: Model,
    pub mesh: Mesh,
    pub control: Control,
}


/**
 * App state
 */
#[derive(Clone, Serialize, Deserialize)]
pub struct App {
    pub state: AnyState,
    pub tasks: Tasks,
    pub config: Configuration,
    pub version: String,
}




// ============================================================================
impl From<State<hydro_euler::euler_2d::Conserved>> for AnyState {
    fn from(state: State<hydro_euler::euler_2d::Conserved>) -> Self {
        Self::Newtonian(state)
    }
}

impl From<State<hydro_srhd::srhd_2d::Conserved>> for AnyState {
    fn from(state: State<hydro_srhd::srhd_2d::Conserved>) -> Self {
        Self::Relativistic(state)
    }
}

impl From<NewtonianHydro> for AnyHydro {
    fn from(hydro: NewtonianHydro) -> Self {
        Self::Newtonian(hydro)
    }
}

impl From<RelativisticHydro> for AnyHydro {
    fn from(hydro: RelativisticHydro) -> Self {
        Self::Relativistic(hydro)
    }
}

impl AnyHydro {
    pub fn validate(&self) -> anyhow::Result<()> {
        match self {
            AnyHydro::Newtonian(hydro) => hydro.validate(),
            AnyHydro::Relativistic(hydro) => hydro.validate(),
        }        
    }
}

impl Control {
    pub fn validate(&self) -> anyhow::Result<()> {
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
impl InitialModel for Model {
    fn validate(&self) -> anyhow::Result<()> {
        match self {
            Model::JetInCloud(m)   => m.validate(),
            Model::HaloKilonova(m) => m.validate(),
        }
    }

    fn primitive_at(&self, coordinate: (f64, f64), time: f64) -> AnyPrimitive {
        match self {
            Model::JetInCloud(m)   => m.primitive_at(coordinate, time),
            Model::HaloKilonova(m) => m.primitive_at(coordinate, time),
        }
    }

    fn scalar_at(&self, coordinate: (f64, f64), time: f64) -> f64 {
        match self {
            Model::JetInCloud(m)   => m.scalar_at(coordinate, time),
            Model::HaloKilonova(m) => m.scalar_at(coordinate, time),
        }
    }
}




// ============================================================================
impl Configuration {
    pub fn package<H>(hydro: &H, model: &Model, mesh: &Mesh, control: &Control) -> Self
    where
        H: Hydrodynamics,
        AnyHydro: From<H> {
        Configuration {
            hydro: AnyHydro::from(hydro.clone()),
            model: model.clone(),
            mesh: mesh.clone(),
            control: control.clone(),
        }
    }

    pub fn validate(&self) -> anyhow::Result<()> {
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
    pub fn validate(self) -> anyhow::Result<Self> {
        self.config.validate()?;
        Ok(self)
    }

    /**
     * Construct a new App instance from a user configuration.
     */
    pub fn from_config(mut config: Configuration) -> Result<Self, Error> {
        for extra_config_str in std::env::args().skip_while(|s| !s.contains('=')) {
            if extra_config_str.ends_with(".yaml") {
                config.patch_from_reader(File::open(extra_config_str)?)?
            } else {
                config.patch_from_key_val(&extra_config_str)?
            }
        }

        let geometry = config.mesh.grid_blocks_geometry(config.control.start_time);
        let state = match &config.hydro {
            AnyHydro::Newtonian(hydro) => {
                AnyState::from(State::from_model(&config.model, hydro, &geometry, config.control.start_time))
            },
            AnyHydro::Relativistic(hydro) => {
                AnyState::from(State::from_model(&config.model, hydro, &geometry, config.control.start_time))
            },
        };
        let tasks = Tasks::new(config.control.start_time);
        Ok(Self{state, tasks, config, version: VERSION_AND_BUILD.to_string()})
    }

    /**
     * Construct a new App instance from a file: may be a config.yaml or a
     * chkpt.0000.cbor.
     */
    pub fn from_file(filename: &str) -> Result<Self, Error> {
        match Path::new(&filename).extension().and_then(OsStr::to_str) {
            Some("yaml") => Self::from_config(serde_yaml::from_str(&read_to_string(filename)?)?),
            Some("cbor") => Ok(io::read_cbor(filename)?),
            _ => Err(Error::UnknownInputType(filename.to_string())),
        }
    }

    /**
     * Construct a new App instance from a preset (hard-coded) configuration
     * name, or otherwise an input file if no matching preset is found.
     */
    pub fn from_preset_or_file(input: &str) -> Result<Self, Error> {
        match input {
            "jet_in_cloud" => Self::from_config(serde_yaml::from_str(std::include_str!("../setups/jet_in_cloud.yaml"))?),
            _ => Self::from_file(input),
        }
    }

    /**
     * Construct a new App instance from references to the member variables.
     */
    pub fn package<C, H>(state: &State<C>, tasks: &mut Tasks, hydro: &H, model: &Model, mesh: &Mesh, control: &Control) -> Self
    where
        H: Hydrodynamics<Conserved = C>,
        C: Conserved,
        AnyState: From<State<C>>,
        AnyHydro: From<H>
    {
        Self {
            state: AnyState::from(state.clone()),
            tasks: tasks.clone(),
            config: Configuration::package(hydro, model, mesh, control),
            version: VERSION_AND_BUILD.to_string(),
        }
    }
}
