use crate::dataset;
use clap::{Args, Parser, Subcommand, ValueEnum};
use serde::Serialize;
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

/// Rebar CLI dataset arguments
#[derive(Parser, Debug)]
pub struct DatasetArgs {
    #[clap(subcommand)]
    pub command: DatasetCommand,
}

/// Rebar CLI dataset command (download, list)
#[derive(Subcommand, Debug)]
#[clap(verbatim_doc_comment)]
pub enum DatasetCommand {
    /// List available Rebar datasets.
    List(DatasetListArgs),

    /// Download available Rebar datasets.
    Download(DatasetDownloadArgs),
}

// -----------------------------------------------------------------------------
// Dataset Get

/// Download available Rebar datasets.
#[derive(Parser, Debug)]
#[clap(verbatim_doc_comment)]
pub struct DatasetDownloadArgs {
    /// Dataset name.
    #[clap(short = 'r', long, required = true)]
    pub name: dataset::Name,

    /// Dataset tag.
    #[clap(short = 't', long)]
    #[clap(default_value_t=dataset::Tag::default())]
    pub tag: dataset::Tag,

    /// Output directory.
    ///
    /// If the directory does not exist, it will be created.
    #[clap(short = 'o', long, required = true)]
    pub output_dir: PathBuf,
}

// -----------------------------------------------------------------------------
// Dataset List

#[derive(Parser, Debug)]
#[clap(verbatim_doc_comment)]
#[group(id = "outputs", required = true, multiple = false)]
/// Download available Rebar dataset.
pub struct DatasetListArgs {
    /// List Rebar datasets with this name.
    #[clap(short = 'n', long)]
    pub name: Option<String>,

    /// List Rebar datasets with this tag.
    #[clap(short = 't', long)]
    #[clap(default_value = "latest")]
    pub tag: String,
}

// -----------------------------------------------------------------------------
// Run
// -----------------------------------------------------------------------------

/// Run Rebar on input alignment or dataset population.
#[derive(Clone, Debug, Parser)]
#[clap(verbatim_doc_comment)]
pub struct RunArgs {
    #[command(flatten)]
    pub input: RunInput,

    /// Dataset directory.
    #[clap(short = 'd', long, required = true)]
    pub dataset_dir: PathBuf,

    /// Number of bases to mask at the 5' and 3' ends.
    #[arg(short = 'm', long, default_value_t = 200)]
    pub mask: usize,

    /// Maximum number of parents.
    #[arg(short = 'p', long, default_value_t = 2)]
    pub max_parents: usize,

    /// Maximum number of search iterations to find parents.
    #[arg(short = 'i', long, default_value_t = 3)]
    pub max_iter: usize,

    /// Minimum number of consecutive bases in a parental region.
    #[arg(short = 'c', long, default_value_t = 3)]
    pub min_consecutive: usize,

    /// Minimum length of a parental region.
    #[arg(short = 'l', long, default_value_t = 500)]
    pub min_length: usize,

    /// Minimum number of substitutions in a parental region.
    #[arg(short = 's', long, default_value_t = 1)]
    pub min_subs: usize,

    /// Output directory.
    ///
    /// If the directory does not exist, it will be created.
    #[clap(short = 'o', long, required = true)]
    pub output_dir: PathBuf,
}

#[derive(Args, Clone, Debug)]
#[group(required = true, multiple = true)]
pub struct RunInput {
    /// Input fasta alignment
    #[arg(long)]
    pub populations: Option<String>,

    /// Input dataset population
    #[arg(long)]
    pub alignment: Option<PathBuf>,
}

// -----------------------------------------------------------------------------
// Plot
// -----------------------------------------------------------------------------

/// Plot Rebar output from input barcodes tsv or barcodes directory.
#[derive(Clone, Debug, Parser)]
#[clap(verbatim_doc_comment)]
pub struct PlotArgs {
    #[command(flatten)]
    pub barcodes: PlotInput,

    /// Input linelist from rebar run.
    #[clap(short = 'i', long, required = true)]
    pub linelist: PathBuf,

    /// Output directory.
    ///
    /// If the directory does not exist, it will be created.
    #[clap(short = 'o', long, required = true)]
    pub output_dir: PathBuf,
}

#[derive(Args, Clone, Debug)]
#[group(required = true, multiple = true)]
pub struct PlotInput {
    /// Input barcodes tsv
    #[arg(short = 'f', long)]
    pub barcodes_file: Option<PathBuf>,

    /// Input barcodes directory
    #[arg(short = 'd', long)]
    pub barcodes_dir: Option<PathBuf>,
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
