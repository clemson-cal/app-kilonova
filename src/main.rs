use kilonova::*;
use app::{
    AgnosticHydro,
    AgnosticState,
    App,
    Configuration,
    Control,
    Model,
};
use mesh::{
    Mesh,
};
use products::{
    Products,
};
use state::{
    State,
};
use traits::{
    Conserved,
    Hydrodynamics,
};
use tasks::{
    Tasks,
};




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
        let products = Products::from_state(state, hydro, &config);
        io::write_cbor(&products, &filename, control.snappy_compression)?;
    }

    if tasks.write_checkpoint.next_time <= state.time {
        tasks.write_checkpoint.advance(control.checkpoint_interval);
        let filename = format!("{}/chkpt.{:04}.{}", outdir, tasks.write_checkpoint.count - 1, extension);
        let app = App::package(state, tasks, hydro, model, mesh, control);
        io::write_cbor(&app, &filename, control.snappy_compression)?;
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
    let outdir = io::parent_directory(&input);

    println!();
    println!("\t{}", app::DESCRIPTION);
    println!("\t{}", app::VERSION_AND_BUILD);
    println!();
    println!("\tinput file ........ {}", input);
    println!("\toutput drectory ... {}", outdir);

    let App{state, tasks, config, ..} = App::from_preset_or_file(&input)?.validate()?;

    for line in serde_yaml::to_string(&config)?.split("\n").skip(1) {
        println!("\t{}", line);
    }
    println!();

    let Configuration{hydro, model, mesh, control} = config;

    match (state, hydro) {
        (AgnosticState::Newtonian(state), AgnosticHydro::Newtonian(hydro)) => {
            run(state, tasks, hydro, model, mesh, control, outdir)
        },
        (AgnosticState::Relativistic(state), AgnosticHydro::Relativistic(hydro)) => {
            run(state, tasks, hydro, model, mesh, control, outdir)
        },
        _ => unreachable!(),
    }
}
