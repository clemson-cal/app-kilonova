use std::collections::HashMap;
use clap::Clap;
use io_logical::verified;




// ============================================================================
pub static DESCRIPTION: &str = env!("CARGO_PKG_DESCRIPTION");
pub static VERSION_AND_BUILD: &str = git_version::git_version!(prefix=concat!("v", env!("CARGO_PKG_VERSION"), " "));




// ============================================================================
#[derive(Clap)]
#[clap(version=VERSION_AND_BUILD, author=clap::crate_authors!(", "))]
pub struct App
{
    #[clap(about="Model parameters")]
    model_parameters: Vec<String>,

    #[clap(short, long, about="Restart file or directory [use latest checkpoint if directory]")]
    restart: Option<String>,

    #[clap(short, long, about="Output directory [default: data/ or restart directory]")]
    outdir: Option<String>,

    #[clap(long, default_value="1", about="Number of iterations between side effects")]
    fold: usize,

    #[clap(long, default_value="1", about="Number of worker threads to use")]
    threads: usize,
}




// ============================================================================
impl App
{
    fn restart_file(&self) -> anyhow::Result<Option<verified::File>>
    {
        if let Some(restart) = self.restart.clone() {
            Ok(Some(verified::file_or_most_recent_matching_in_directory(restart, "chkpt.????.h5")?))
        } else {
            Ok(None)
        }
    }

    fn output_rundir_child(&self, filename: &str) -> anyhow::Result<Option<verified::File>>
    {
        if let Ok(file) = self.output_directory()?.existing_child(filename) {
            Ok(Some(file))
        } else {
            Ok(None)
        }
    }

    fn output_directory(&self) -> anyhow::Result<verified::Directory>
    {
        if let Some(outdir) = self.outdir.clone() {
            Ok(verified::Directory::require(outdir)?)
        } else if let Some(restart) = &self.restart_file()? {
            Ok(restart.parent())
        } else {
            Ok(verified::Directory::require("data".into())?)
        }
    }

    fn restart_model_parameters(&self) -> anyhow::Result<HashMap<String, kind_config::Value>>
    {
        if let Some(file) = self.restart_file()? {
        	Ok(kind_config::io::read_from_hdf5(&hdf5::File::open(file.as_str())?.group("model")?)?)
        } else {
            Ok(HashMap::new())
        }
    }

    fn compute_units(&self, num_blocks: usize) -> usize
    {
        num_cpus::get_physical().min(self.threads)
    }
}
