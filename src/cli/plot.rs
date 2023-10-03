use clap::Parser;
use std::path::PathBuf;

/// Plot Rebar output from input barcodes tsv or barcodes directory.
#[derive(Clone, Debug, Parser)]
#[clap(verbatim_doc_comment)]
pub struct Args {
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
