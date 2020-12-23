#![allow(unused)]
use clap::Clap;




mod app;
mod physics;
mod scheme;
mod state;
mod tasks;
mod traits;




fn main() {
	let app = app::App::parse();
}
