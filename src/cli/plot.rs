use clap::Parser;
use std::path::PathBuf;

/// Plot recombination from 'run' output.
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

    /// Directory to download fonts to.
    #[clap(short = 'f', long, default_value = PathBuf::from(".cache/fonts").into_os_string())]
    pub font_cache: PathBuf,
}
