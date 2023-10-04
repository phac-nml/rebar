use crate::cli::run;
use crate::dataset::attributes::Summary;
use crate::dataset::Dataset;
use crate::phylogeny::Phylogeny;
use crate::sequence::{read_reference, Sequence, Substitution};
use crate::utils::table::Table;
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use log::{debug, info, warn};
use std::collections::BTreeMap;
use std::path::Path;
use std::str::FromStr;

// ----------------------------------------------------------------------------
// Dataset
// ----------------------------------------------------------------------------

/// Load dataset.
pub fn dataset(args: &run::Args) -> Result<Dataset, Report> {
    info!("Loading dataset: {:?}", args.dataset_dir);

    let mut dataset = Dataset::new();

    // Summary
    let summary_path = args.dataset_dir.join("summary.json");
    let summary = Summary::read(&summary_path)?;
    dataset.name = summary.name;
    dataset.tag = summary.tag;

    // Reference
    let reference_path = args.dataset_dir.join("reference.fasta");
    dataset.reference = read_reference(&reference_path, &args.mask)?;

    // Populations and Mutations
    let populations_path = args.dataset_dir.join("populations.fasta");
    (dataset.populations, dataset.mutations) =
        parse_populations(&populations_path, &reference_path, &args.mask)?;

    // Edge Cases
    let edge_cases_path = args.dataset_dir.join("edge_cases.json");
    let multiple = true;
    dataset.edge_cases = run::Args::read(&edge_cases_path, multiple)?.unwrap_right();

    // ------------------------------------------------------------------------
    // Load Phylogeny (Optional)

    let phylogeny_path = args.dataset_dir.join("phylogeny.json");
    dataset.phylogeny = if phylogeny_path.exists() {
        Phylogeny::read(&phylogeny_path)?
    } else {
        warn!("No phylogeny was found.");
        Phylogeny::new()
    };

    // --------------------------------------------------------------------
    // Diagnostic Mutations (Optional)

    let diagnostic_path = args.dataset_dir.join("diagnostic_mutations.tsv");
    dataset.diagnostic = BTreeMap::new();

    if diagnostic_path.exists() {
        let diagnostic_table = Table::read(&diagnostic_path)?;
        let mut_col_i = diagnostic_table.header_position("mutation")?;
        let pop_col_i = diagnostic_table.header_position("population")?;

        for row in diagnostic_table.rows {
            let mutation = Substitution::from_str(&row[mut_col_i])?;
            let population = row[pop_col_i].clone();
            dataset.diagnostic.insert(mutation, population);
        }
    } else {
        warn!("No diagnostic mutations were found.");
    }

    // --------------------------------------------------------------------
    // Done

    Ok(dataset)
}

// ----------------------------------------------------------------------------
// Parse Populations
// ----------------------------------------------------------------------------

#[allow(clippy::type_complexity)]
pub fn parse_populations(
    populations_path: &Path,
    reference_path: &Path,
    mask: &Vec<usize>,
) -> Result<
    (
        BTreeMap<String, Sequence>,
        BTreeMap<Substitution, Vec<String>>,
    ),
    Report,
> {
    // read in populations from fasta
    let populations_reader = fasta::Reader::from_file(populations_path)
        .map_err(|e| eyre!(e))
        .wrap_err("Failed to read file: {populations_path:?}")?;

    // read in reference from fasta
    let reference = read_reference(reference_path, mask)?;

    let mut populations = BTreeMap::new();
    let mut mutations = BTreeMap::new();

    for result in populations_reader.records() {
        let record = result?;
        let sequence = Sequence::from_record(record, Some(&reference), mask)?;
        populations.insert(sequence.id.clone(), sequence.clone());
        for sub in sequence.substitutions {
            mutations
                .entry(sub)
                .or_insert(Vec::new())
                .push(sequence.id.clone());
        }
    }

    Ok((populations, mutations))
}

pub fn get_diagnostic_mutations(
    mutations: &BTreeMap<Substitution, Vec<String>>,
    phylogeny: &Phylogeny,
) -> Result<Table, Report> {
    let mut table = Table::new();
    table.headers = vec!["mutation", "population", "include_descendants"]
        .into_iter()
        .map(String::from)
        .collect_vec();

    // configure progress bar style
    let progress_bar_style = ProgressStyle::with_template(
        "{bar:40} {pos}/{len} ({percent}%) | Mutations / Second: {per_sec} | Elapsed: {elapsed_precise} | ETA: {eta_precise}"
    ).unwrap();
    let progress_bar = ProgressBar::new(mutations.len() as u64);
    progress_bar.set_style(progress_bar_style);

    for (mutation, populations) in mutations {
        progress_bar.inc(1);
        let mut population = populations[0].clone();
        debug!("{mutation:?}");
        debug!("\tpopulations: {}", populations.len());
        let mut is_diagnostic = false;
        let mut include_descendants = false;

        // A mutation found in only 1 population is obviously diagnostic
        if populations.len() == 1 {
            is_diagnostic = true;
        } else if !phylogeny.is_empty() {
            let common_ancestor = phylogeny.get_common_ancestor(populations)?;
            if populations.contains(&common_ancestor) {
                is_diagnostic = true;
                population = common_ancestor;
                include_descendants = true;
            }
        }

        if is_diagnostic {
            let row = vec![
                mutation.to_string(),
                population,
                include_descendants.to_string(),
            ];
            table.rows.push(row);
        }
    }

    progress_bar.inc(1);

    Ok(table)
}
