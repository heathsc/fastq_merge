#[macro_use]
extern crate log;
#[macro_use]
extern crate clap;

use std::io;

pub mod config;
pub mod cli;
pub mod stats;
mod process;

fn main() -> io::Result<()> {
	let config = cli::handle_cli()?;
	debug!("Starting processing");
	process::process(&config)	
}
