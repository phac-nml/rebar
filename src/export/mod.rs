use crate::dataset::{Dataset, SearchResult};
use crate::recombination::Recombination;
use crate::sequence::Sequence;
use color_eyre::eyre::{Report, Result};
use csv;
use itertools::Itertools;
use log::info;
use std::collections::BTreeMap;
use std::path::Path;

pub fn collect_recombinants(linelist: &Vec<Vec<String>>) -> Result<(), Report> {
    // a unique recombinant is defined by the same population(?), parent, breakpoints
    let mut uniq_recombinants: BTreeMap<String, Vec<String>> = BTreeMap::new();

    // read in headers
    let headers = linelist.iter().next().unwrap();
    info!("headers: {headers:?}");

    let strain_col = headers.iter().position(|h| h == "strain").unwrap();
    let recombinant_col = headers.iter().position(|h| h == "recombinant").unwrap();
    let parents_col = headers.iter().position(|h| h == "parents").unwrap();
    let breakpoints_col = headers.iter().position(|h| h == "breakpoints").unwrap();

    for row in linelist {
        info!("row: {row:?}");
        let strain = &row[strain_col];
        let recombinant = &row[recombinant_col];
        let parents = &row[parents_col];
        let breakpoints = &row[breakpoints_col];

        let uniq_key = format!("{}_{}_{}", recombinant, parents, breakpoints);
        uniq_recombinants
            .entry(uniq_key.to_string())
            .or_insert(Vec::new())
            .push(strain.to_string());
    }

    Ok(())
}

// pub fn write_barcodes(
//     output_dir: &Path,
//     sequences: &Vec<Sequence>,
//     recombinations: &BTreeMap<String, Recombination>,
// ) -> Result<(), Report> {
//     Ok(())
// }

pub fn write_linelist(
    output_path: &Path,
    sequences: &[Sequence],
    best_matches: &BTreeMap<String, SearchResult>,
    recombinations: &BTreeMap<String, Recombination>,
    dataset: &Dataset,
) -> Result<Vec<Vec<String>>, Report> {
    // setup table (vector of vectors)
    let mut table = Vec::new();

    // table headers, convert from &str to String
    let headers = vec![
        "strain",
        "population",
        "recombinant",
        "parents",
        "breakpoints",
        "regions",
        "genome_length",
        "dataset_name",
        "dataset_tag",
    ]
    .iter()
    .map(|h| h.to_string())
    .collect_vec();

    table.push(headers.clone());

    for seq in sequences.iter() {
        if !best_matches.contains_key(&seq.id) || !recombinations.contains_key(&seq.id) {
            continue;
        }

        let best_match = &best_matches[&seq.id];
        let recombination = &recombinations[&seq.id];

        let mut row = vec!["".to_string(); headers.len()];

        // strain
        row[0] = seq.id.to_string();
        // population
        row[1] = best_match.consensus_population.to_string();
        // recombinant
        let recombinant = &best_match.recombinant;
        if let Some(recombinant) = recombinant {
            row[2] = recombinant.to_string();
        }
        // parents
        row[3] = recombination.parents.join(",").to_string();
        // breakpoints
        row[4] = recombination.breakpoints.iter().join(",").to_string();
        // regions
        row[5] = recombination.regions.values().join(",").to_string();
        // genome_length
        row[6] = seq.genome_length.to_string();
        // dataset name
        row[7] = dataset.name.to_string();
        // dataset tag
        row[8] = dataset.tag.to_string();

        // add row to table
        table.push(row);
    }

    // write table to file
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output_path)?;

    for row in &table {
        writer.write_record(row)?;
    }

    Ok(table)
}
