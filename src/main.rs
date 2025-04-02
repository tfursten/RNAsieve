mod cli;
mod build;
mod filter;
// mod utils;

use clap::Parser;
use cli::{Cli, Commands};

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::BuildIndex(args) => {
            if let Err(e) = build::build_index(args) {
                eprintln!("Error: {}", e);
            }
        },
        Commands::Filter(args) => {
            if let Err(e) = filter::filter(args) {
                eprintln!("Error: {}", e);
            }
        },
}}

