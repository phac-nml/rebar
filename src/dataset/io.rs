use crate::cli;
use crate::dataset::attributes::{Name, Summary, SummaryExportFormat, Tag};
use crate::dataset::edge_cases::{EdgeCaseExportFormat, EdgeCaseImportFormat};
use crate::dataset::{edge_cases, sarscov2, Dataset};
use crate::phylogeny::{Phylogeny, PhylogenyExportFormat, PhylogenyImportFormat};
use crate::sequence::{read_reference, Sequence, Substitution};
use crate::utils;
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result};
use indicatif::{style::ProgressStyle, ProgressBar};
use log::{debug, info, warn};
use std::collections::BTreeMap;
use std::fs::{create_dir_all, File};
use std::path::Path;
use std::str::FromStr;

// ----------------------------------------------------------------------------
// Dataset Download
// ----------------------------------------------------------------------------

/// Download a remote dataset
pub async fn download_dataset(args: &cli::DatasetDownloadArgs) -> Result<(), Report> {
    // Create the output directory if it doesn't exist
    if !args.output_dir.exists() {
        create_dir_all(&args.output_dir)?;
        info!("Creating output directory: {:?}", args.output_dir);
    }

    // --------------------------------------------------------------------
    // Download Reference

    let url = match args.name {
        Name::SarsCov2 => sarscov2::REFERENCE_URL.to_string(),
        _ => {
            return Err(eyre!(
                "Reference download for {} is not implemented.",
                args.name
            ))
        }
    };
    let ext = utils::path_to_ext(Path::new(&url))?;
    let mut decompress = false;
    if ext != "fasta" && ext != "fa" {
        decompress = true;
    }
    let reference_path = args.output_dir.join("reference.fasta");
    info!("Downloading reference: {} to {:?}", url, reference_path);
    utils::download_file(&url, &reference_path, decompress).await?;

    // --------------------------------------------------------------------
    // Create Annotations

    let annotations_path = args.output_dir.join("annotations.tsv");
    info!("Downloading annotations to {:?}", annotations_path);

    let annotations = match args.name {
        Name::SarsCov2 => sarscov2::create_annotations()?,
        _ => {
            return Err(eyre!(
                "Annotations for {} dataset is not implemented.",
                args.name
            ))
        }
    };
    let annotations_table = annotations.to_table()?;
    annotations_table.write(&annotations_path)?;

    // --------------------------------------------------------------------
    // Download Populations

    let url = match args.name {
        Name::SarsCov2 => sarscov2::POPULATIONS_URL.to_string(),
        _ => {
            return Err(eyre!(
                "Population download for {} is not implemented.",
                args.name
            ))
        }
    };
    let ext = utils::path_to_ext(Path::new(&url))?;
    let mut decompress = false;
    if ext != "fasta" && ext != "fa" {
        decompress = true;
    }
    let populations_path = args.output_dir.join("populations.fasta");
    info!("Downloading populations: {} to {:?}", url, populations_path);
    utils::download_file(&url, &populations_path, decompress).await?;

    // --------------------------------------------------------------------
    // Create Phylogeny

    let phylogeny_path = args.output_dir.join("phylogeny.json");
    info!("Creating phylogeny: {:?}", phylogeny_path);

    let mut phylogeny = Phylogeny::new();
    phylogeny.build_graph(&args.name, &args.output_dir).await?;
    phylogeny.export(&args.output_dir, PhylogenyExportFormat::Dot)?;
    phylogeny.export(&args.output_dir, PhylogenyExportFormat::Json)?;

    // --------------------------------------------------------------------
    // Identify Diagnostic Mutations
    //
    // This is painfully slow, need to rethink!

    if args.diagnostic {
        let diagnostic_path = args.output_dir.join("diagnostic_mutations.tsv");
        info!("Identifying diagnostic mutations: {:?}", diagnostic_path);

        let mask = 0;
        let (_populations, mutations) =
            parse_populations(&populations_path, &reference_path, mask)?;
        let diagnostic_table = get_diagnostic_mutations(&mutations, &phylogeny)?;
        diagnostic_table.write(&diagnostic_path)?;
    }

    // --------------------------------------------------------------------
    // Create Edge Cases

    let edge_cases_path = args.output_dir.join("edge_cases.json");
    info!("Creating edge cases: {:?}", edge_cases_path);

    let edge_cases = match args.name {
        Name::SarsCov2 => sarscov2::create_edge_cases()?,
        _ => {
            return Err(eyre!(
                "Edge cases for {} dataset is not implemented.",
                args.name
            ))
        }
    };
    edge_cases::export(&edge_cases, &args.output_dir, EdgeCaseExportFormat::Json)?;

    // --------------------------------------------------------------------
    // Create Summary

    let output_path = args.output_dir.join("summary.json");
    info!("Creating info summary: {:?}", output_path);

    let summary = Summary {
        name: args.name,
        tag: args.tag.clone(),
    };
    summary.export(&args.output_dir, SummaryExportFormat::Json)?;

    // --------------------------------------------------------------------
    // Finish

    info!("Done.");
    Ok(())
}

// ----------------------------------------------------------------------------
// Dataset Load
// ----------------------------------------------------------------------------

/// Load a local dataset
pub fn load_dataset(args: &cli::RunArgs) -> Result<Dataset, Report> {
    // Initialize a new container for our dataset pieces
    let mut dataset = Dataset::new();

    // ------------------------------------------------------------------------
    // Load Reference (Required)

    let reference_path = args.dataset_dir.join("reference.fasta");
    info!("Loading reference: {reference_path:?}");
    dataset.reference = read_reference(&reference_path, args.mask)?;

    // ------------------------------------------------------------------------
    // Load Populations and Parse Mutations (Required)

    let populations_path = args.dataset_dir.join("populations.fasta");
    info!("Loading populations: {populations_path:?}");
    (dataset.populations, dataset.mutations) =
        parse_populations(&populations_path, &reference_path, args.mask)?;

    // ------------------------------------------------------------------------
    // Load Attributes/Summary (optional)

    let summary_path = args.dataset_dir.join("summary.json");
    if summary_path.exists() {
        info!("Loading summary: {:?}", summary_path);
        let reader = File::open(summary_path)?;
        let summary: Summary = serde_json::from_reader(&reader)?;
        dataset.name = summary.name;
        dataset.tag = summary.tag;
    } else {
        warn!("Optional summary file was not found: {summary_path:?}");
        dataset.tag = Tag::Unknown;
        dataset.name = Name::Unknown;
    }

    // ------------------------------------------------------------------------
    // Load Phylogeny (Optional)

    let phylogeny_path = args.dataset_dir.join("phylogeny.json");
    dataset.phylogeny = Phylogeny::new();
    if phylogeny_path.exists() {
        info!("Loading phylogeny: {phylogeny_path:?}");
        dataset.phylogeny =
            Phylogeny::import(&args.dataset_dir, PhylogenyImportFormat::Json)?;
    } else {
        warn!("Optional phylogeny file was not found: {phylogeny_path:?}");
    }

    // --------------------------------------------------------------------
    // Load Diagnostic Mutations (Optional)

    let diagnostic_path = args.dataset_dir.join("diagnostic_mutations.tsv");
    dataset.diagnostic = BTreeMap::new();
    if diagnostic_path.exists() {
        info!("Loading diagnostic mutations: {diagnostic_path:?}");
        let diagnostic_table = utils::read_table(&diagnostic_path)?;
        let mut_col_i = diagnostic_table.header_position("mutation")?;
        let pop_col_i = diagnostic_table.header_position("population")?;

        for row in diagnostic_table.rows {
            let mutation = Substitution::from_str(&row[mut_col_i])?;
            let population = row[pop_col_i].clone();
            dataset.diagnostic.insert(mutation, population);
        }
    } else {
        warn!("Optional diagnostic mutations file was not found: {diagnostic_path:?}");
    }

    // --------------------------------------------------------------------
    // Load/Create Edge Cases (Optional)

    let edge_cases_path = args.dataset_dir.join("edge_cases.json");
    dataset.edge_cases = Vec::new();

    if edge_cases_path.exists() {
        info!("Loading edge cases: {edge_cases_path:?}");
        dataset.edge_cases =
            edge_cases::import(&args.dataset_dir, EdgeCaseImportFormat::Json)?
    } else {
        warn!("Optional edge cases file was not found: {edge_cases_path:?}");
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
) -> Result<utils::table::Table, Report> {
    let mut table = utils::table::Table::new();
    table.headers = vec!["mutation", "population", "include_descendants"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

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
