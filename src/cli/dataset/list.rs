use crate::dataset::attributes::Name;
use clap::Parser;

// -----------------------------------------------------------------------------
// Dataset List

/// List datasets.
#[derive(Parser, Debug)]
#[clap(verbatim_doc_comment)]
pub struct Args {
    /// Dataset name.
    #[clap(short = 'n', long)]
    pub name: Option<Name>,
}
