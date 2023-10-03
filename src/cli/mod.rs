pub mod dataset;
pub mod run;
pub mod plot;

use clap::{Parser, Subcommand, ValueEnum};
use serde::Serialize;
use std::default::Default;

// -----------------------------------------------------------------------------
// CLI Entry Point
// -----------------------------------------------------------------------------

/// Rebar command-line interface (CLI)
#[derive(Parser, Debug)]
#[clap(name = "rebar", trailing_var_arg = true)]
#[clap(author, version)]
#[clap(verbatim_doc_comment)]
pub struct Cli {
    #[clap(subcommand)]
    // rebar command (dataset, run, help)
    pub command: Command,

    /// Control output verbosity level.
    #[clap(short = 'v', long)]
    #[clap(value_enum, default_value_t = Verbosity::default())]
    #[clap(hide_possible_values = false)]
    #[clap(global = true)]
    pub verbosity: Verbosity,
}

/// Rebar CLI commands (dataset, run, plot).
#[derive(Subcommand, Debug)]
#[clap(verbatim_doc_comment)]
pub enum Command {
    Dataset(Box<dataset::Args>),
    Run(Box<run::Args>),
    Plot(Box<plot::Args>),
}

// -----------------------------------------------------------------------------
// Verbosity
// -----------------------------------------------------------------------------

#[derive(Clone, Debug, Default, Serialize, ValueEnum)]
pub enum Verbosity {
    #[default]
    Info,
    Warn,
    Debug,
    Error,
}

impl std::fmt::Display for Verbosity {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Convert to lowercase for RUST_LOG env var compatibility
        let lowercase = format!("{:?}", self).to_lowercase();
        write!(f, "{lowercase}")
    }
}