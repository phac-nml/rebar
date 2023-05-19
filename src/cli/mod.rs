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

    /// Make output more quiet or more verbose
    #[arg(short, long)]
    pub verbosity: Option<String>,
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
    },

    /// Run rebar on input alignment.
    Run {
        /// Dataset directory
        #[arg(short, long, required = true)]
        dataset_dir: PathBuf,

        /// Input fasta alignment
        #[arg(short, long, required = true)]
        alignment: PathBuf,

        /// Mask 5' and 3' ends of alignment
        #[arg(short, long, default_value_t = 200)]
        mask: usize,
    },
}
