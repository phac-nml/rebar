use crate::cli::run;
use crate::dataset::attributes::{Name, Summary, Tag};
use crate::dataset::Dataset;
use crate::phylogeny::Phylogeny;
use crate::sequence::{read_reference, Sequence, Substitution};
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use log::{info, warn};
use std::collections::BTreeMap;
use std::path::Path;

// ----------------------------------------------------------------------------
// Dataset
// ----------------------------------------------------------------------------

/// Load dataset.
pub fn dataset(dataset_dir: &Path, mask: &Vec<usize>) -> Result<Dataset, Report> {
    info!("Loading dataset: {:?}", dataset_dir);

    let mut dataset = Dataset::new();

    // ------------------------------------------------------------------------
    // Mandatory

    // Reference
    let reference_path = dataset_dir.join("reference.fasta");
    dataset.reference = read_reference(&reference_path, mask)?;

    // Populations and Mutations
    let populations_path = dataset_dir.join("populations.fasta");
    (dataset.populations, dataset.mutations) =
        parse_populations(&populations_path, &reference_path, mask)?;

    // ------------------------------------------------------------------------
    // Optional

    // Summary
    let summary_path = dataset_dir.join("summary.json");
    if summary_path.exists() {
        let summary = Summary::read(&summary_path)?;
        dataset.name = summary.name;
        dataset.tag = summary.tag;
    } else {
        warn!("No summary was found: {summary_path:?}");
        dataset.name = Name::Custom;
        dataset.tag = Tag::Custom;
    }

    // Edge Cases
    let edge_cases_path = dataset_dir.join("edge_cases.json");
    dataset.edge_cases = if edge_cases_path.exists() {
        let multiple = true;
        run::Args::read(&edge_cases_path, multiple)?.unwrap_right()
    } else {
        warn!("No edge cases were found: {edge_cases_path:?}");
        Vec::new()
    };

    // Phylogeny
    let phylogeny_path = dataset_dir.join("phylogeny.json");
    dataset.phylogeny = if phylogeny_path.exists() {
        Phylogeny::read(&phylogeny_path)?
    } else {
        warn!("No phylogeny was found: {phylogeny_path:?}");
        Phylogeny::new()
    };

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
        .wrap_err(format!("Failed to read file: {populations_path:?}"))?;

    // read in reference from fasta
    let reference = read_reference(reference_path, mask)?;

    let mut populations = BTreeMap::new();
    let mut mutations = BTreeMap::new();

    for result in populations_reader.records() {
        let record = result?;
        let sequence = Sequence::from_record(record, Some(&reference), mask)?;
        populations.insert(sequence.id.clone(), sequence.clone());

        for sub in sequence.substitutions {
            mutations.entry(sub).or_insert(Vec::new()).push(sequence.id.clone());
        }
    }

    Ok((populations, mutations))
}
