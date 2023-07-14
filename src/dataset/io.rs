use crate::annotation::Annotations;
use crate::cli;
use crate::dataset::edge_cases::{EdgeCaseExportFormat, EdgeCaseImportFormat};
use crate::dataset::{attributes, constants, edge_cases, Dataset};
use crate::phylogeny::{Phylogeny, PhylogenyExportFormat, PhylogenyImportFormat};
use crate::sequence::{Sequence, Substitution};
use crate::utils;
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result};
use log::info;
use std::collections::BTreeMap;
use std::fs::{create_dir_all, File};
use std::path::Path;
use std::str::FromStr;

// ----------------------------------------------------------------------------
// Dataset Download
// ----------------------------------------------------------------------------

/// Download a remote dataset
pub async fn download(
    name: &attributes::Name,
    tag: &attributes::Tag,
    output_dir: &Path,
) -> Result<(), Report> {
    if !output_dir.exists() {
        create_dir_all(output_dir)?;
        info!("Creating output directory: {:?}", output_dir);
    }

    // --------------------------------------------------------------------
    // Download Reference

    let url = match name {
        attributes::Name::SarsCov2 => constants::SARSCOV2_REFERENCE_URL.to_string(),
        _ => {
            return Err(eyre!(
                "Downloading the {name} dataset is not implemented yet."
            ))
        }
    };
    let ext = Path::new(&url).extension().unwrap().to_str().unwrap();
    let mut decompress = false;
    if ext != "fasta" && ext != "fa" {
        decompress = true;
    }
    let reference_path = output_dir.join("reference.fasta");
    info!("Downloading reference: {} to {:?}", url, reference_path);
    utils::download_file(&url, &reference_path, decompress).await?;

    // --------------------------------------------------------------------
    // Create Annotations

    let annotations_path = output_dir.join("annotations.tsv");
    info!("Downloading annotations to {:?}", annotations_path);
    let annotations = Annotations::from_name(name)?;
    let annotations_table = annotations.to_table()?;
    annotations_table.write(&annotations_path, Some('\t'))?;

    // --------------------------------------------------------------------
    // Download Populations

    let url = match name {
        attributes::Name::SarsCov2 => constants::SARSCOV2_POPULATIONS_URL.to_string(),
        _ => {
            return Err(eyre!(
                "Downloading the {name} dataset is not implemented yet."
            ))
        }
    };
    let ext = Path::new(&url).extension().unwrap().to_str().unwrap();
    let mut decompress = false;
    if ext != "fasta" && ext != "fa" {
        decompress = true;
    }
    let populations_path = output_dir.join("populations.fasta");
    info!("Downloading populations: {} to {:?}", url, populations_path);
    utils::download_file(&url, &populations_path, decompress).await?;

    // --------------------------------------------------------------------
    // Create Phylogeny

    let phylogeny_path = output_dir.join("phylogeny.dot");
    info!("Creating phylogeny: {:?}", phylogeny_path);

    let mut phylogeny = Phylogeny::new();
    phylogeny.build_graph(name, output_dir).await?;
    phylogeny.export(output_dir, PhylogenyExportFormat::Dot)?;
    phylogeny.export(output_dir, PhylogenyExportFormat::Json)?;

    // --------------------------------------------------------------------
    // Identify Diagnostic Mutations
    //
    // This is painfully slow, need to rethink!

    // let diagnostic_path = output_dir.join("diagnostic_mutations.tsv");
    // info!("Identifying diagnostic mutations: {:?}", diagnostic_path);

    // let mask = 0;
    // let (_populations, mutations) = parse_populations(&populations_path, &reference_path, mask)?;
    // let diagnostic_table = identify_diagnostic_mutations(&mutations, &phylogeny)?;
    // diagnostic_table.write(&diagnostic_path, Some('\t'))?;

    // --------------------------------------------------------------------
    // Create Edge Cases

    let edge_cases_path = output_dir.join("edge_cases.yaml");
    info!("Creating edge cases: {:?}", edge_cases_path);

    let edge_cases = edge_cases::from_dataset_name(name)?;
    edge_cases::export(&edge_cases, output_dir, EdgeCaseExportFormat::Json)?;

    // --------------------------------------------------------------------
    // Create Summary

    let output_path = output_dir.join("summary.yaml");
    info!("Creating info summary: {:?}", output_path);

    let summary = attributes::Summary {
        name: *name,
        tag: tag.clone(),
    };
    summary.export(output_dir, attributes::SummaryExportFormat::Json)?;

    // --------------------------------------------------------------------
    // Finish

    info!("Done.");
    Ok(())
}

// ----------------------------------------------------------------------------
// Dataset Load
// ----------------------------------------------------------------------------

/// Load a local dataset
pub fn load(args: &cli::RunArgs) -> Result<Dataset, Report> {
    // Initialize a new container for our dataset pieces
    let mut dataset = Dataset::new();

    // ------------------------------------------------------------------------
    // Load Reference (Required)

    let reference_path = args.dataset_dir.join("reference.fasta");
    info!("Loading reference: {:?}", reference_path);
    // read in reference from fasta
    let reference_reader = fasta::Reader::from_file(reference_path.clone())
        .expect("Unable to load reference");
    // parse just the first record (ie. next() )
    let reference = reference_reader.records().next().unwrap().unwrap();
    // convert to a bio Sequence object
    dataset.reference = Sequence::from_record(reference, None, args.mask)?;

    // ------------------------------------------------------------------------
    // Load Populations and Parse Mutations (Required)

    let populations_path = args.dataset_dir.join("populations.fasta");
    (dataset.populations, dataset.mutations) =
        parse_populations(&populations_path, &reference_path, args.mask)?;

    // ------------------------------------------------------------------------
    // Load Attributes/Summary (optional)

    let summary_path = args.dataset_dir.join("summary.json");
    if summary_path.exists() {
        info!("Loading summary: {:?}", summary_path);
        let reader = File::open(summary_path)?;
        let summary: attributes::Summary = serde_json::from_reader(&reader)?;
        dataset.name = summary.name;
        dataset.tag = summary.tag;
    } else {
        dataset.tag = attributes::Tag::Unknown;
        dataset.name = attributes::Name::Unknown;
    }

    // ------------------------------------------------------------------------
    // Load Phylogeny (Optional)

    let phylogeny_path = args.dataset_dir.join("phylogeny.json");
    dataset.phylogeny = Phylogeny::new();
    if phylogeny_path.exists() {
        info!("Loading phylogeny: {:?}", phylogeny_path);
        dataset.phylogeny =
            Phylogeny::import(&args.dataset_dir, PhylogenyImportFormat::Json)?;
    }

    // --------------------------------------------------------------------
    // Load Diagnostic Mutations (Optional)

    let diagnostic_path = args.dataset_dir.join("diagnostic_mutations.tsv");
    dataset.diagnostic = BTreeMap::new();
    if diagnostic_path.exists() {
        info!("Loading diagnostic mutations: {diagnostic_path:?}");
        let diagnostic_table = utils::Table::from_tsv(&diagnostic_path)?;
        let mut_col_i = diagnostic_table.header_position("mutation")?;
        let pop_col_i = diagnostic_table.header_position("population")?;

        for row in diagnostic_table.rows {
            let mutation = Substitution::from_str(&row[mut_col_i])?;
            let population = row[pop_col_i].clone();
            dataset.diagnostic.insert(mutation, population);
        }
    }

    // --------------------------------------------------------------------
    // Load/Create Edge Cases (Optional)

    let edge_cases_path = args.dataset_dir.join("edge_cases.yaml");
    dataset.edge_cases = if edge_cases_path.exists() {
        info!("Loading edge cases: {edge_cases_path:?}");
        edge_cases::import(&args.dataset_dir, EdgeCaseImportFormat::Json)?
    } else {
        Vec::new()
    };

    // --------------------------------------------------------------------
    // Done

    Ok(dataset)
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

#[allow(clippy::type_complexity)]
pub fn parse_populations(
    populations_path: &Path,
    reference_path: &Path,
    mask: usize,
) -> Result<
    (
        BTreeMap<String, Sequence>,
        BTreeMap<Substitution, Vec<String>>,
    ),
    Report,
> {
    // read in reference from fasta
    let populations_reader =
        fasta::Reader::from_file(populations_path).expect("Unable to load populations");

    // read in reference from fasta
    let reference_reader =
        fasta::Reader::from_file(reference_path).expect("Unable to load reference");
    // parse just the first record (ie. next() )
    let reference = reference_reader.records().next().unwrap().unwrap();
    // convert to a bio Sequence object
    let reference = Sequence::from_record(reference, None, mask)?;

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

pub fn identify_diagnostic_mutations(
    mutations: &BTreeMap<Substitution, Vec<String>>,
    phylogeny: &Phylogeny,
) -> Result<utils::Table, Report> {
    let mut table = utils::Table::new();
    table.headers = vec!["mutation", "population", "include_descendants"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    for (mutation, populations) in mutations {
        info!("mutation: {mutation:?}");

        let mut population = populations[0].clone();

        let mut is_diagnostic = false;
        let mut include_descendants = false;

        // A mutation found in only 1 population is obviously diagnostic
        if populations.len() == 1 {
            is_diagnostic = true;
        }
        // Otherwise with the phylogeny, a diagnostic mutation is monophyletic
        // but may not include all descendants
        // Note: I'm not 100% sure of this logic at the moment
        else if !phylogeny.is_empty() {
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

    Ok(table)
}
