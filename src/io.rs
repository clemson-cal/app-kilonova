use std::collections::HashMap;
use hdf5::Group;
use hdf5::File;
use io_logical::verified;
use io_logical::nicer_hdf5;
use io_logical::nicer_hdf5::{H5Read, H5Write};
use crate::app::VERSION_AND_BUILD;
use crate::tasks::Tasks;
use crate::state::State;
use crate::traits::Conserved;




// ============================================================================
impl nicer_hdf5::H5Read for Tasks {
    fn read(group: &Group, name: &str) -> hdf5::Result<Self> {
        nicer_hdf5::read_as_keyed_vec(group, name)
    }    
}

impl nicer_hdf5::H5Write for Tasks {
    fn write(&self, group: &Group, name: &str) -> hdf5::Result<()> {
        nicer_hdf5::write_as_keyed_vec(self.clone(), group, name)
    }
}




// ============================================================================
fn write_state<C: Conserved>(group: &Group, state: &State<C>) -> hdf5::Result<()> {
    let state_group = group.create_group("state")?;
    let solution_group = state_group.create_group("solution")?;

    for (index, state) in &state.solution {
        let block_group = solution_group.create_group(&index.to_string())?;
        state.conserved.write(&block_group, "conserved")?;
        state.scalar.write(&block_group, "scalar")?;
    }
    state.time.write(&state_group, "time")?;
    state.iteration.write(&state_group, "iteration")?;
    Ok(())
}

fn write_build(group: &Group) -> hdf5::Result<()> {
    use hdf5::types::VarLenAscii;

    group.new_dataset::<VarLenAscii>()
        .create("version", ())?
        .write_scalar(&VarLenAscii::from_ascii(VERSION_AND_BUILD).unwrap())?;

    Ok(())
}




// ============================================================================
fn write_tasks(group: &Group, tasks: &Tasks) -> hdf5::Result<()> {
    tasks.write(group, "tasks")
}

fn read_tasks(file: &verified::File) -> hdf5::Result<Tasks> {
    let file = File::open(file.as_str())?;
    Tasks::read(&file, "tasks")
}

fn write_model(group: &Group, model: &kind_config::Form) -> hdf5::Result<()> {
    kind_config::io::write_to_hdf5(&group.create_group("model")?, &model.value_map())
}

fn read_model(file: &verified::File) -> hdf5::Result<HashMap::<String, kind_config::Value>> {
    kind_config::io::read_from_hdf5(&File::open(file.as_str())?.group("model")?)
}




// ============================================================================
pub fn write_checkpoint<C: Conserved>(
    filename: &str,
    state: &State<C>,
    model: &kind_config::Form,
    tasks: &Tasks) -> hdf5::Result<()> {
    let file = File::create(filename)?;

    write_state(&file, &state)?;
    write_tasks(&file, &tasks)?;
    write_model(&file, &model)?;
    write_build(&file)?;

    Ok(())
}
