use crate::utils::table::Table;
use color_eyre::eyre::{Report, Result};
use itertools::Itertools;

/// Create Toy1 genome annotations.
pub fn build() -> Result<Table, Report> {
    let mut table = Table::new();

    let headers = vec!["gene", "abbreviation", "start", "end"];
    let rows = vec![
        vec!["Gene1", "g1", "1", "3"],
        vec!["Gene2", "g2", "12", "20"],
    ];

    // Convert values to String
    table.headers = headers.into_iter().map(String::from).collect_vec();
    table.rows = rows
        .into_iter()
        .map(|row| row.into_iter().map(String::from).collect_vec())
        .collect_vec();

    Ok(table)
}
