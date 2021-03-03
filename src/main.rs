use kilonova::*;
use app::{
    AnyHydro,
    AnyModel,
    AnyState,
    App,
    Configuration,
    Control,
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
    InitialModel,
};
use tasks::{
    Tasks,
};




// ============================================================================
fn side_effects<C, M, H>(state: &State<C>, tasks: &mut Tasks, hydro: &H, model: &M, mesh: &Mesh, control: &Control)
    -> anyhow::Result<()>
where
    H: Hydrodynamics<Conserved = C>,
    M: InitialModel,
    C: Conserved,
    AnyHydro: From<H>,
    AnyModel: From<M>,
    AnyState: From<State<C>>,
{
    if tasks.iteration_message.next_time <= state.time {
        let time = tasks.iteration_message.advance(0.0);
        let mzps = 1e-6 * state.total_zones() as f64 / time * control.fold as f64;
        if tasks.iteration_message.count_this_run > 1 {
            println!("[{:05}] t={:.5} blocks={} Mzps={:.2})", state.iteration, state.time, state.solution.len(), mzps);
        }
    }

    if let Some(products_interval) = control.products_interval {
        if tasks.write_products.next_time <= state.time {
            tasks.write_products.advance(products_interval);
            let filename = format!("{}/prods.{:04}.cbor", control.output_directory, tasks.write_products.count - 1);
            let config = Configuration::package(hydro, model, mesh, control);
            let products = Products::try_from_state(state, hydro, &config)?;
            std::fs::create_dir_all(&control.output_directory)?;
            io::write_cbor(&products, &filename)?;
        }
    }

    if tasks.write_checkpoint.next_time <= state.time {
        tasks.write_checkpoint.advance(control.checkpoint_interval);
        let filename = format!("{}/chkpt.{:04}.cbor", control.output_directory, tasks.write_checkpoint.count - 1);
        let app = App::package(state, tasks, hydro, model, mesh, control);
        std::fs::create_dir_all(&control.output_directory)?;
        io::write_cbor(&app, &filename)?;
    }

    Ok(())
}




// ============================================================================
fn run<C, M, H>(mut state: State<C>, mut tasks: Tasks, hydro: H, model: M, mesh: Mesh, control: Control)
    -> anyhow::Result<()>
where
    H: Hydrodynamics<Conserved = C>,
    M: InitialModel,
    C: Conserved,
    AnyHydro: From<H>,
    AnyModel: From<M>,
    AnyState: From<State<C>>,
{
    let mut block_geometry = mesh.grid_blocks_geometry(state.time);
    let runtime = tokio::runtime::Builder::new_multi_thread()
        .worker_threads(control.num_threads())
        .build()?;

    while state.time < control.final_time {
        side_effects(&state, &mut tasks, &hydro, &model, &mesh, &control)?;
        state = scheme::advance(state, &hydro, &model, &mesh, &mut block_geometry, &runtime, control.fold)?;
    }

    side_effects(&state, &mut tasks, &hydro, &model, &mesh, &control)?;

    Ok(())
}




// ============================================================================
fn main() -> anyhow::Result<()> {

    println!();
    println!("{}", app::DESCRIPTION);
    println!("{}", app::VERSION_AND_BUILD);
    println!();

    match std::env::args().nth(1) {
        None => {
            println!("usage: kilonova <input.yaml|chkpt.cbor|preset> [opts.yaml|group.key=value] [...]");
            println!();
            println!("These are the preset model setups:");
            println!();
            for (key, _) in App::presets() {
                println!("  {}", key);
            }
            println!();
            println!("To run any of these presets, run e.g. `kilonova jet_in_star`.");
            Ok(())
        }
        Some(input) => {
            let overrides = std::env::args().skip(2).collect();
            let App{state, tasks, config, ..} = App::from_preset_or_file(&input, overrides)?.validate()?;

            for line in serde_yaml::to_string(&config)?.split("\n").skip(1) {
                println!("{}", line);
            }
            println!();

            let Configuration{hydro, model, mesh, control} = config;

            println!("worker threads ...... {}", control.num_threads());
            println!("compute cores ....... {}", num_cpus::get());
            println!();

            match (state, hydro) {
                (AnyState::Newtonian(state), AnyHydro::Newtonian(hydro)) => {
                    run(state, tasks, hydro, model, mesh, control)
                },
                (AnyState::Relativistic(state), AnyHydro::Relativistic(hydro)) => {
                    run(state, tasks, hydro, model, mesh, control)
                },
                _ => unreachable!(),
            }
        }
    }
}
