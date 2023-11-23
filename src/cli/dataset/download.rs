use crate::dataset::attributes::{Name, Tag};
use clap::Parser;
use std::path::PathBuf;

/// Download dataset.
#[derive(Parser, Debug)]
#[clap(verbatim_doc_comment)]
pub struct Args {
    /// Dataset name.
    #[clap(short = 'r', long, required = true)]
    pub name: Name,

    /// Dataset tag.
    ///
    /// can be 'latest' or a date (YYYY-MM-DD)
    #[clap(short = 't', long, required = true)]
    pub tag: Tag,

    /// Output directory.
    ///
    /// If the directory does not exist, it will be created.
    #[clap(short = 'o', long, required = true)]
    pub output_dir: PathBuf,

    /// Download dataset from a summary.json snapshot.
    #[clap(short = 's', long)]
    pub summary: Option<PathBuf>,
}
