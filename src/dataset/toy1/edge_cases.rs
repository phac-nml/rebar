use crate::cli::run;
use color_eyre::eyre::{Report, Result};

/// Create default Toy1 recombinant edge cases.
pub fn default() -> Result<Vec<run::Args>, Report> {
    let edge_cases: Vec<run::Args> = Vec::new();
    Ok(edge_cases)
}
