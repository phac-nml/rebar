use crate::dataset::{Dataset, SearchResult};
use crate::recombination;
use crate::utils;
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

// ----------------------------------------------------------------------------
// Validate

#[derive(Clone, Debug, Deserialize, Serialize)]
pub enum Validate {
    // population status
    CorrectPopulation,
    IncorrectPopulation,
    // // recombinant status
    // PositiveRecombinant,
    // NegativeRecombinant,
    // FalsePositiveRecombinant,
    // FalseNegativeRecombinant,
    // // parents status
    // CorrectParents,
    // IncorrectParents,
}

impl FromStr for Validate {
    type Err = Report;

    fn from_str(result: &str) -> Result<Self, Report> {
        let validate = match result {
            "correct_population" => Validate::CorrectPopulation,
            "incorrect_population" => Validate::IncorrectPopulation,
            _ => return Err(eyre!("Unknown validation result: {result}")),
        };

        Ok(validate)
    }
}

impl fmt::Display for Validate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let validate = match self {
            Validate::CorrectPopulation => "correct_population",
            Validate::IncorrectPopulation => "incorrect_population",
        };

        write!(f, "{}", validate)
    }
}

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
        "population",
        "recombinant",
        "parents",
        "breakpoints",
        "validate",
        "unique_key",
        "regions",
        "private",
        "genome_length",
        "dataset_name",
        "dataset_tag",
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
        if best_match.recombinant.is_some() {
            let recombinant = best_match.recombinant.clone().unwrap();
            row[table.header_position("recombinant")?] = recombinant;
        }

        // parents
        let parents = recombination.parents.join(",").to_string();
        row[table.header_position("parents")?] = parents;

        // breakpoints
        let breakpoints = recombination.breakpoints.iter().join(",").to_string();
        row[table.header_position("breakpoints")?] = breakpoints;

        // unique_key
        let unique_key = recombination.unique_key.to_string();
        row[table.header_position("unique_key")?] = unique_key;

        // regions
        let regions = recombination.regions.values().join(",").to_string();
        row[table.header_position("regions")?] = regions;

        // private mutations
        let private = best_match.private.iter().join(",").to_string();
        row[table.header_position("private")?] = private;

        // genome_length
        let genome_length = recombination.genome_length.to_string();
        row[table.header_position("genome_length")?] = genome_length;

        // dataset name
        row[table.header_position("dataset_name")?] = dataset.name.to_string();

        // dataset tag
        row[table.header_position("dataset_tag")?] = dataset.tag.to_string();

        // --------------------------------------------------------------------
        // Validate

        let mut validate = Vec::new();

        let strain = strain.replace("query_", "");

        if dataset.populations.contains_key(&strain) {
            // Population status
            if strain == population {
                validate.push(Validate::CorrectPopulation);
            } else {
                validate.push(Validate::IncorrectPopulation);
            }
            // Recombinant status
            // Parent status
        }

        row[table.header_position("validate")?] = validate.iter().join(",").to_string();

        table.rows.push(row);
    }

    Ok(table)
}
