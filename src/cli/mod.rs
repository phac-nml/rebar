pub mod verbosity;

use clap::{Parser, Subcommand};
use serde::Serialize;
use std::path::PathBuf;

#[derive(Debug, Parser, Serialize)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
#[clap(verbatim_doc_comment)]
/// Recombination barcode detector.
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand, Serialize)]
pub enum Command {
    /// Download and list available datasets.
    Dataset {
        /// Dataset name
        #[arg(short, long, required = true)]
        name: String,

        /// Dataset tag
        #[arg(short, long, required = true)]
        tag: String,

        // /// Dataset reference accession
        // #[arg(short, long, required = true)]
        // reference: String,
        /// Dataset output directory
        #[arg(short, long, required = true)]
        output_dir: PathBuf,

        /// Make output more quiet or more verbose
        #[arg(short, long, default_value_t=String::from("info"))]
        verbosity: String,
    },

    /// Run rebar on input alignment.
    Run {
        /// Dataset directory
        #[arg(short = 'd', long, required = true)]
        dataset_dir: PathBuf,

        /// Input fasta alignment
        #[arg(short = 'a', long, required = true)]
        alignment: PathBuf,

        /// Mask 5' and 3' ends of alignment
        #[arg(short = 'm', long, default_value_t = 200)]
        mask: usize,

        /// Maximum number of parents
        #[arg(short = 'p', long, default_value_t = 2)]
        max_parents: usize,

        /// Make output more quiet or more verbose
        #[arg(short, long, default_value_t=String::from("info"))]
        verbosity: String,
    },
}
