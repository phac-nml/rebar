use clap::Parser;
use crate::dataset::attributes::{Name, Tag};
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
    /// can be 'latest' or a date in RFC 3339 format (2023-09-01T12:00:00T)
    #[clap(short = 't', long)]
    #[clap(default_value_t=Tag::default())]
    pub tag: Tag,

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
