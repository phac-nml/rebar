use clap::Parser;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Run on input alignment or dataset population.
#[derive(Clone, Debug, Deserialize, Parser, Serialize)]
#[clap(verbatim_doc_comment)]
pub struct Args {
    /// Dataset directory.
    #[clap(short = 'd', long, required = true)]
    pub dataset_dir: PathBuf,

    /// Simulate recombination between these parents.
    ///
    /// 5' -> 3'
    #[arg(long, value_delimiter = ',', required = true)]
    pub parents: Vec<String>,

    /// Specify the breakpoints.
    ///
    /// If not provided, will be random.
    #[arg(long, value_delimiter = ',')]
    pub breakpoints: Option<Vec<usize>>,

    /// Output directory.
    ///
    /// If the directory does not exist, it will be created.
    #[clap(short = 'o', long, required = true)]
    pub output_dir: PathBuf,
}

impl Default for Args {
    fn default() -> Self {
        Self::new()
    }
}

impl Args {
    pub fn new() -> Self {
        Args {
            breakpoints: None,
            dataset_dir: PathBuf::new(),
            output_dir: PathBuf::new(),
            parents: Vec::new(),
        }
    }
}
