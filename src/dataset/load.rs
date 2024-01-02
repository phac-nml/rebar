use crate::cli::run;
use crate::dataset::attributes::{Name, Summary, Tag};
use crate::dataset::Dataset;
use crate::phylogeny::Phylogeny;
use crate::sequence::{read_reference, Sequence, Substitution};
use color_eyre::eyre::{Report, Result};
use log::{info, warn};
use noodles::fasta;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
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
    let mut reader =
        File::open(populations_path).map(BufReader::new).map(fasta::Reader::new)?;

    // read in reference from fasta
    let reference = read_reference(reference_path, mask)?;

    let mut populations = BTreeMap::new();
    let mut mutations = BTreeMap::new();

    reader.records().try_for_each(|result| {
        let record = result?;
        let sequence = Sequence::from_record(record, Some(&reference), mask)?;
        populations.insert(sequence.id.clone(), sequence.clone());

        sequence.substitutions.into_iter().for_each(|sub| {
            mutations.entry(sub).or_insert(Vec::new()).push(sequence.id.clone());
        });
        Ok::<(), Report>(())
    })?;

    Ok((populations, mutations))
}
