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
    #[clap(short = 'r', long, required = true)]
    pub run_dir: PathBuf,

    /// A single barcodes TSV file to plot.
    #[clap(short = 'b', long)]
    pub barcodes_file: Option<PathBuf>,

    /// Output directory for plots.
    ///
    /// Otherwise will default to 'plots/' under the --run-dir
    #[clap(short = 'o', long)]
    pub output_dir: Option<PathBuf>,

    /// Draw all coordinates in barcodes file.
    ///
    /// By default, rebar wil only draw coordinates where the parent populations have different bases.
    /// As a result, private mutations at these coordinates will not be visualized.
    /// This flag will forcibly draw all coordinates in the barcodes file, but be warned, depending on
    /// the genome size and number of samples, this may cause a crash.
    #[clap(short = 'p', long)]
    pub all_coords: bool,
}

impl Default for Args {
    fn default() -> Self {
        Self::new()
    }
}

impl Args {
    pub fn new() -> Self {
        Args {
            dataset_dir: PathBuf::new(),
            run_dir: PathBuf::new(),
            barcodes_file: None,
            output_dir: None,
            all_coords: false,
        }
    }
}
