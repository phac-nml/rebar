use crate::dataset;
use clap::{Args, Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};
use std::default::Default;
use std::path::PathBuf;

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
    Dataset(Box<DatasetArgs>),
    Run(Box<RunArgs>),
    Plot(Box<PlotArgs>),
}

// -----------------------------------------------------------------------------
// Dataset
// -----------------------------------------------------------------------------

/// Dataset arguments
#[derive(Parser, Debug)]
pub struct DatasetArgs {
    #[clap(subcommand)]
    pub command: DatasetCommand,
}

/// Dataset commands (download, list)
#[derive(Subcommand, Debug)]
#[clap(verbatim_doc_comment)]
pub enum DatasetCommand {
    /// List datasets.
    List(DatasetListArgs),

    /// Download dataset.
    Download(DatasetDownloadArgs),
}

// -----------------------------------------------------------------------------
// Dataset Get

/// Download dataset.
#[derive(Parser, Debug)]
#[clap(verbatim_doc_comment)]
pub struct DatasetDownloadArgs {
    /// Dataset name.
    #[clap(short = 'r', long, required = true)]
    pub name: dataset::attributes::Name,

    /// Dataset tag.
    ///
    /// can be 'latest' or a date in RFC 3339 format (2023-09-01T12:00:00T)
    #[clap(short = 't', long)]
    #[clap(default_value_t=dataset::attributes::Tag::default())]
    pub tag: dataset::attributes::Tag,

    /// Calculate diagnostic mutations (slow).
    #[clap(short = 'd', long)]
    #[clap(default_value_t = false)]
    pub diagnostic: bool,

    /// Output directory.
    ///
    /// If the directory does not exist, it will be created.
    #[clap(short = 'o', long, required = true)]
    pub output_dir: PathBuf,

    /// Download dataset from a summary.json snapshot.
    #[clap(short = 's', long)]
    pub summary: Option<PathBuf>,
}

// -----------------------------------------------------------------------------
// Dataset List

/// List datasets.
#[derive(Parser, Debug)]
#[clap(verbatim_doc_comment)]
pub struct DatasetListArgs {
    /// Dataset name.
    #[clap(short = 'n', long)]
    pub name: Option<dataset::attributes::Name>,
}

// -----------------------------------------------------------------------------
// Run
// -----------------------------------------------------------------------------

/// Run on input alignment or dataset population.
#[derive(Clone, Debug, Deserialize, Parser, Serialize)]
#[clap(verbatim_doc_comment)]
pub struct RunArgs {
    #[command(flatten)]
    pub input: RunInput,

    /// Dataset directory.
    #[clap(short = 'd', long, required = true)]
    pub dataset_dir: PathBuf,

    /// Number of bases to mask at the 5' and 3' ends.
    #[arg(short = 'm', long, default_value_t = RunArgs::default().mask)]
    pub mask: usize,

    /// Maximum number of parents.
    #[arg(short = 'p', long, default_value_t = RunArgs::default().max_parents)]
    pub max_parents: usize,

    /// Maximum number of search iterations to find each parent.
    #[arg(short = 'i', long, default_value_t = RunArgs::default().max_iter)]
    pub max_iter: usize,

    /// Minimum number of consecutive bases in a parental region.
    #[arg(short = 'c', long, default_value_t = RunArgs::default().min_consecutive)]
    pub min_consecutive: usize,

    /// Minimum length of a parental region.
    #[arg(short = 'l', long, default_value_t = RunArgs::default().min_length)]
    pub min_length: usize,

    /// Minimum number of substitutions in a parental region.
    #[arg(short = 's', long, default_value_t = RunArgs::default().min_subs)]
    pub min_subs: usize,

    /// Output directory.
    ///
    /// If the directory does not exist, it will be created.
    #[clap(short = 'o', long, required = true)]
    pub output_dir: PathBuf,

    /// Number of CPU threads to use.
    #[clap(short = 't', long, default_value_t = RunArgs::default().threads)]
    pub threads: usize,

    /// Run an unbiased search, which does not use information about designated recombinant parents.
    #[arg(short = 'u', long, default_value_t = RunArgs::default().unbiased)]
    pub unbiased: bool,
}

impl Default for RunArgs {
    fn default() -> Self {
        RunArgs {
            input: RunInput::default(),
            dataset_dir: PathBuf::new(),
            mask: 200,
            max_parents: 2,
            max_iter: 3,
            min_consecutive: 3,
            min_length: 500,
            min_subs: 1,
            output_dir: PathBuf::new(),
            threads: 1,
            unbiased: false,
        }
    }
}

impl RunArgs {
    pub fn new() -> Self {
        RunArgs {
            input: RunInput::default(),
            dataset_dir: PathBuf::new(),
            mask: 0,
            max_parents: 0,
            max_iter: 0,
            min_consecutive: 0,
            min_length: 0,
            min_subs: 0,
            output_dir: PathBuf::new(),
            threads: 0,
            unbiased: false,
        }
    }
}

#[derive(Args, Clone, Debug, Deserialize, Serialize)]
#[group(required = true, multiple = true)]
pub struct RunInput {
    /// Input fasta alignment
    #[arg(long)]
    pub populations: Option<String>,

    /// Input dataset population
    #[arg(long)]
    pub alignment: Option<PathBuf>,
}

impl Default for RunInput {
    fn default() -> Self {
        Self::new()
    }
}

impl RunInput {
    pub fn new() -> Self {
        RunInput {
            populations: None,
            alignment: None,
        }
    }
}

// -----------------------------------------------------------------------------
// Plot
// -----------------------------------------------------------------------------

/// Plot Rebar output from input barcodes tsv or barcodes directory.
#[derive(Clone, Debug, Parser)]
#[clap(verbatim_doc_comment)]
pub struct PlotArgs {
    /// rebar dataset directory.
    #[clap(short = 'd', long, required = true)]
    pub dataset_dir: PathBuf,

    /// Output directory from rebar run.
    ///
    /// Will plot all TSV files under barcodes/
    #[clap(short = 'o', long, required = true)]
    pub output_dir: PathBuf,

    /// A single barcodes TSV file to plot.
    #[clap(short = 'b', long)]
    pub barcodes_file: Option<PathBuf>,

    /// Output directory for plots.
    ///
    /// Otherwise will default to 'plots' under the --output-dir
    #[clap(short = 'p', long)]
    pub plot_dir: Option<PathBuf>,
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
