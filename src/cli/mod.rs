pub mod log;

use clap::{Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct Cli {

    #[command(subcommand)]
    command: Option<Commands>,

    // /// Name of the person to greet
    // #[arg(short, long)]
    // name: String,

    // /// Number of times to greet
    // #[arg(short, long, default_value_t = 1)]
    // count: u8,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Download and list available datasets.
    Dataset { name: Option<String> },
    /// Run rebar on querry sequences.
    Run { name: Option<String> },
}
