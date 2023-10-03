pub mod download;
pub mod list;

use clap::{Parser, Subcommand};

/// Dataset arguments
#[derive(Parser, Debug)]
pub struct Args {
    #[clap(subcommand)]
    pub command: Command,
}

/// Dataset commands (download, list)
#[derive(Subcommand, Debug)]
#[clap(verbatim_doc_comment)]
pub enum Command {
    /// List datasets.
    List(list::Args),

    /// Download dataset.
    Download(download::Args),
}
