use crate::dataset::{Dataset, SearchResult};
use crate::recombination::{validate, Recombination};
use crate::utils;
use color_eyre::eyre::{Report, Result};
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

        let subs_by_origin = recombination.get_substitution_origins(best_match)?;
        let mut origins = Vec::new();
        // origin order: primary parent, secondary parent, recombinant, private
        if recombination.recombinant.is_some() {
            origins.extend(recombination.parents.clone());
        } else {
            origins.push(best_match.consensus_population.clone());
        }
        origins.push("private".to_string());

        // string format
        let substitutions = origins
            .iter()
            .filter_map(|o| {
                let subs = subs_by_origin.get(o).cloned().unwrap_or_default();
                let subs_format = format!("{}|{o}", subs.iter().join(","));
                (!subs.is_empty()).then_some(subs_format)
            })
            .join(";");
        row[table.header_position("substitutions")?] = substitutions;

        table.rows.push(row);
    }

    Ok(table)
}
