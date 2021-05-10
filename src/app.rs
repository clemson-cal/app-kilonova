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
    HaloKilonova,
    JetInCloud,
    JetInStar,
    WindShock,
    KineticBomb,
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
#[derive(Clone, Serialize, Deserialize, derive_more::From)]
#[serde(deny_unknown_fields, rename_all = "snake_case")]
pub enum AnyModel {
    HaloKilonova(HaloKilonova),
    JetInCloud(JetInCloud),
    JetInStar(JetInStar),
    WindShock(WindShock),
    KineticBomb(KineticBomb),
}




/**
 * Enum for any of the supported hydrodynamics types
 */
#[derive(Clone, Serialize, Deserialize, derive_more::From)]
#[serde(deny_unknown_fields, rename_all = "snake_case")]
pub enum AnyHydro {
    Newtonian(NewtonianHydro),
    Relativistic(RelativisticHydro),
}




/**
 * Enum for the solution state of any of the supported hydrodynamics types
 */
#[derive(Clone, Serialize, Deserialize, derive_more::From)]
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

    /// The simulation start time. This is not necessarily t=0, because
    /// model setups may have a time-dependent background solution.
    pub start_time: f64,

    /// The simulation end time.
    pub final_time: f64,

    /// The time between writing checkpoint  files.
    pub checkpoint_interval: f64,

    /// The time between writing products files. If omitted or nil, defaults
    /// to no products output. This option should be considered deprecated.
    /// Write checkpoints and then convert them to products files in
    /// post-processing if needed.
    pub products_interval: Option<f64>,

    /// The number of iterations between performing side-effects
    pub fold: usize,

    /// Number of worker threads on the Tokio runtime. If omitted or nil,
    /// defaults to 2x the number of physical cores.
    pub num_threads: Option<usize>,

    /// Deprecated
    #[serde(default)]
    pub snappy_compression: bool,

    /// The directory where data file will be output. If omitted or nil,
    /// defaults to a the current directory.
    #[serde(default = "Control::default_output_directory")]
    pub output_directory: String,
}

impl Control {
    pub fn num_threads(&self) -> usize {
        match self.num_threads {
            Some(n) => n,
            None => num_cpus::get() * 2,
        }
    }
    fn default_output_directory() -> String {
        ".".into()
    }
}




/**
 * User configuration
 */
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Configuration {
    pub hydro: AnyHydro,
    pub model: AnyModel,
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
        if self.num_threads() == 0 || self.num_threads() >= 1024 {
            anyhow::bail!("num_threads must be > 0 and < 1024")
        }
        if self.checkpoint_interval < 0.0 {
            anyhow::bail!("checkpoint_interval <= 0.0")
        }
        if self.products_interval.unwrap_or(0.0) < 0.0 {
            anyhow::bail!("products_interval <= 0.0")
        }
        Ok(())
    }
}




// ============================================================================
impl InitialModel for AnyModel {

    fn validate(&self) -> anyhow::Result<()> {
        match self {
            AnyModel::HaloKilonova(m) => m.validate(),
            AnyModel::JetInCloud(m)   => m.validate(),
            AnyModel::JetInStar(m)    => m.validate(),
            AnyModel::WindShock(m)    => m.validate(),
            AnyModel::KineticBomb(m) => m.validate(),
        }
    }

    fn primitive_at(&self, coordinate: (f64, f64), time: f64) -> AnyPrimitive {
        match self {
            AnyModel::HaloKilonova(m) => m.primitive_at(coordinate, time),
            AnyModel::JetInCloud(m)   => m.primitive_at(coordinate, time),
            AnyModel::JetInStar(m)    => m.primitive_at(coordinate, time),
            AnyModel::WindShock(m)    => m.primitive_at(coordinate, time),
            AnyModel::KineticBomb(m)  => m.primitive_at(coordinate, time),
        } 
    }

    fn scalar_at(&self, coordinate: (f64, f64), time: f64) -> f64 {
        match self {
            AnyModel::HaloKilonova(m) => m.scalar_at(coordinate, time),
            AnyModel::JetInCloud(m)   => m.scalar_at(coordinate, time),
            AnyModel::JetInStar(m)    => m.scalar_at(coordinate, time),
            AnyModel::WindShock(m)    => m.scalar_at(coordinate, time),
            AnyModel::KineticBomb(m)  => m.scalar_at(coordinate, time),
        }
    }
}




// ============================================================================
impl Configuration {
    pub fn package<H, M>(hydro: &H, model: &M, mesh: &Mesh, control: &Control) -> Self
    where
        H: Hydrodynamics,
        M: InitialModel,
        AnyHydro: From<H>,
        AnyModel: From<M>,
    {
        Configuration {
            hydro: hydro.clone().into(),
            model: model.clone().into(),
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

    /**
     * Patch this config struct with inputs from the command line. The inputs
     * can be names of YAML files or key=value pairs.
     */
    pub fn patch_from(&mut self, overrides: Vec<String>) -> Result<(), Error> {
        for extra_config_str in overrides {
            if extra_config_str.ends_with(".yaml") {
                self.patch_from_reader(File::open(extra_config_str)?)?
            } else {
                self.patch_from_key_val(&extra_config_str)?
            }
        }
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
    pub fn from_config(mut config: Configuration, overrides: Vec<String>) -> Result<Self, Error> {

        config.patch_from(overrides)?;

        let geometry = config.mesh.grid_blocks_geometry(config.control.start_time);
        let state = match &config.hydro {
            AnyHydro::Newtonian(hydro) => {
                State::from_model(&config.model, hydro, &geometry, config.control.start_time).into()
            },
            AnyHydro::Relativistic(hydro) => {
                State::from_model(&config.model, hydro, &geometry, config.control.start_time).into()
            },
        };
        let tasks = Tasks::new(config.control.start_time);
        Ok(Self{state, tasks, config, version: VERSION_AND_BUILD.to_string()})
    }

    /**
     * Patch the config struct with inputs from the command line.
     */
    pub fn with_patched_config(mut self, overrides: Vec<String>) -> Result<Self, Error> {
        self.config.patch_from(overrides)?;
        Ok(self)
    }

    /**
     * Construct a new App instance from a file: may be a config.yaml or a
     * chkpt.0000.cbor.
     */
    pub fn from_file(filename: &str, overrides: Vec<String>) -> Result<Self, Error> {
        match Path::new(&filename).extension().and_then(OsStr::to_str) {
            Some("yaml") => Self::from_config(serde_yaml::from_str(&read_to_string(filename)?)?, overrides),
            Some("cbor") => Ok(io::read_cbor::<Self>(filename)?.with_patched_config(overrides)?),
            _ => Err(Error::UnknownInputType(filename.to_string())),
        }
    }

    /**
     * Construct a new App instance from a preset (hard-coded) configuration
     * name, or otherwise an input file if no matching preset is found.
     */
    pub fn from_preset_or_file(input: &str, overrides: Vec<String>) -> Result<Self, Error> {
        for (key, yaml) in Self::presets() {
            if input == key {
                return Ok(Self::from_config(serde_yaml::from_str(yaml)?, overrides)?)
            }
        }
        Self::from_file(input, overrides)
    }

    /**
     * Construct a new App instance from references to the member variables.
     */
    pub fn package<H, M, C>(state: &State<C>, tasks: &Tasks, hydro: &H, model: &M, mesh: &Mesh, control: &Control) -> Self
    where
        H: Hydrodynamics<Conserved = C>,
        M: InitialModel,
        C: Conserved,
        AnyHydro: From<H>,
        AnyModel: From<M>,
        AnyState: From<State<C>>,
    {
        Self {
            state: state.clone().into(),
            tasks: tasks.clone(),
            config: Configuration::package(hydro, model, mesh, control),
            version: VERSION_AND_BUILD.to_string(),
        }
    }

    pub fn presets() -> Vec<(&'static str, &'static str)> {
        vec![
            ("jet_in_cloud", include_str!("setups/jet_in_cloud.yaml")),
            ("jet_in_star", include_str!("setups/jet_in_star.yaml")),
            ("halo_kilonova", include_str!("setups/halo_kilonova.yaml")),
            ("wind_shock", include_str!("setups/wind_shock.yaml")),
            ("kinetic_bomb", include_str!("setups/kinetic_bomb.yaml")),
        ]
    }
}
