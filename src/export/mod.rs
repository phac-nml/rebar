use crate::dataset::{Dataset, SearchResult};
use crate::recombination::{validate, Recombination};
use crate::utils;
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;

// ----------------------------------------------------------------------------
// LineList

pub fn linelist(
    results: &Vec<(SearchResult, Recombination)>,
    dataset: &Dataset,
) -> Result<utils::table::Table, Report> {
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
        "substitutions",
        "genome_length",
        "dataset_name",
        "dataset_tag",
        "cli_version",
    ]
    .into_iter()
    .map(|s| s.to_string())
    .collect_vec();

    // iterate in parallel, checking for same sequence id
    for (best_match, recombination) in results {
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

        // validate, currently requires phylogeny
        if !dataset.phylogeny.is_empty() {
            let validate = validate::validate(dataset, best_match, recombination)?;
            if let Some(validate) = validate {
                row[table.header_position("validate")?] = validate.status.to_string();
                row[table.header_position("validate_details")?] =
                    validate.details.iter().join(";");
            }
        }

        // unique_key
        let unique_key = recombination.unique_key.to_string();
        row[table.header_position("unique_key")?] = unique_key;

        // regions
        let regions = recombination.regions.values().join(",").to_string();
        row[table.header_position("regions")?] = regions;

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

        // --------------------------------------------------------------------
        // Substitutions, annotated by parental origin or private
        // todo!() think about if we want reversions in here or not...

        let mut subs_by_origin = Vec::new();
        // by parent
        if recombination.recombinant.is_some() {
            recombination.parents.iter().for_each(|p| {
                let support = &recombination.support[p];
                if !support.is_empty() {
                    let subs_format = format!("{}|{p}", support.iter().join(","));
                    subs_by_origin.push(subs_format)
                }
            });
        } else {
            let support = &best_match.support[&best_match.consensus_population];

            if !support.is_empty() {
                let subs_format = format!(
                    "{}|{}",
                    support.iter().join(","),
                    &best_match.consensus_population
                );
                subs_by_origin.push(subs_format);
            }
        };

        // private substitutions
        let mut private = if recombination.recombinant.is_some() {
            recombination.private.values().flatten().collect_vec()
        } else {
            best_match.private.iter().collect_vec()
        };
        if !private.is_empty() {
            private.sort();
            let subs_format = format!("{}|private", private.iter().join(","));
            subs_by_origin.push(subs_format);
        }
        row[table.header_position("substitutions")?] = subs_by_origin.iter().join(";");

        table.rows.push(row);
    }

    Ok(table)
}
