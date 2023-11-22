use crate::dataset::{Dataset, SearchResult};
use crate::utils;
use crate::{recombination, recombination::validate};
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;

// ----------------------------------------------------------------------------
// LineList

pub fn linelist(
    recombinations: &Vec<recombination::Recombination>,
    best_matches: &Vec<SearchResult>,
    dataset: &Dataset,
) -> Result<utils::table::Table, Report> {
    // check for same length
    if recombinations.len() != best_matches.len() {
        return Err(eyre!(
            "recombinations and best_matches are different lengths."
        ));
    }

    let mut table = utils::table::Table::new();

    table.headers = vec![
        "strain",
        "validate",
        "validate_details",
        "population",
        "recombinant",
        "parents",
        "breakpoints",
        "edge_case",
        "unique_key",
        "regions",
        "private",
        "diagnostic",
        "genome_length",
        "dataset_name",
        "dataset_tag",
        "cli_version",
        "cli_args",
    ]
    .into_iter()
    .map(|s| s.to_string())
    .collect_vec();

    // iterate in parallel, checking for same sequence id
    for it in recombinations.iter().zip(best_matches.iter()) {
        let (recombination, best_match) = it;

        // check that they're in the correct order
        if recombination.sequence.id != best_match.sequence_id {
            return Err(eyre!(
                "Recombination ID {} is not the same as Best Match ID: {}",
                recombination.sequence.id,
                best_match.sequence_id,
            ));
        }

        // initialize the table row
        let mut row = vec![String::new(); table.headers.len()];

        // strain
        let strain = recombination.sequence.id.to_string();
        row[table.header_position("strain")?] = strain.clone();

        // population
        let population = best_match.consensus_population.to_string();
        row[table.header_position("population")?] = population.clone();

        // recombinant
        if let Some(recombinant) = &recombination.recombinant {
            row[table.header_position("recombinant")?] = recombinant.clone();
        }

        // parents
        let parents = recombination.parents.join(",").to_string();
        row[table.header_position("parents")?] = parents;

        // breakpoints
        let breakpoints = recombination.breakpoints.iter().join(",").to_string();
        row[table.header_position("breakpoints")?] = breakpoints;

        // edge_case
        let edge_case = recombination.edge_case.to_string();
        row[table.header_position("edge_case")?] = edge_case;

        // validate
        let validate = validate::validate(dataset, best_match, recombination)?;
        if let Some(validate) = validate {
            row[table.header_position("validate")?] = validate.status.to_string();
            row[table.header_position("validate_details")?] =
                validate.details.iter().join(";");
        }

        // unique_key
        let unique_key = recombination.unique_key.to_string();
        row[table.header_position("unique_key")?] = unique_key;

        // regions
        let regions = recombination.regions.values().join(",").to_string();
        row[table.header_position("regions")?] = regions;

        // private mutations
        let private = best_match.private.iter().join(",").to_string();
        row[table.header_position("private")?] = private;

        // diagnostic mutations
        // mark this as "NA" if diagnostic mutations were not run for the dataset
        let diagnostic = if dataset.diagnostic.is_empty() {
            "NA".to_string()
        } else {
            best_match.diagnostic.iter().join(",").to_string()
        };
        row[table.header_position("diagnostic")?] = diagnostic;

        // genome_length
        let genome_length = recombination.genome_length.to_string();
        row[table.header_position("genome_length")?] = genome_length;

        // dataset name
        row[table.header_position("dataset_name")?] = dataset.name.to_string();

        // dataset tag
        row[table.header_position("dataset_tag")?] = dataset.tag.to_string();

        // cli version
        row[table.header_position("cli_version")?] =
            env!("CARGO_PKG_VERSION").to_string();

        table.rows.push(row);
    }

    Ok(table)
}
