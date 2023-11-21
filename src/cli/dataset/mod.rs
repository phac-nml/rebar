pub mod download;
pub mod list;

use clap::{Parser, Subcommand};

/// List or download datasets.
#[derive(Parser, Debug)]
#[clap(verbatim_doc_comment)]
pub struct Args {
    #[clap(subcommand)]
    pub command: Command,
}

/// List or download datasets.
#[derive(Subcommand, Debug)]
#[clap(verbatim_doc_comment)]
pub enum Command {
    /// List datasets.
    List(list::Args),

    /// Download dataset.
    Download(download::Args),
}
