use crate::cli;
use crate::dataset::attributes::{
    check_compatibility, Name, Summary, SummaryExportFormat, Tag,
};
use crate::dataset::edge_cases::{EdgeCaseExportFormat, EdgeCaseImportFormat};
use crate::dataset::{edge_cases, sarscov2, Dataset};
use crate::phylogeny::{Phylogeny, PhylogenyExportFormat, PhylogenyImportFormat};
use crate::sequence::{read_reference, Sequence, Substitution};
use crate::utils;
use bio::io::fasta;
use chrono::SecondsFormat;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use indicatif::{style::ProgressStyle, ProgressBar};
use itertools::Itertools;
use log::{debug, info, warn};
use std::collections::BTreeMap;
use std::fs::{create_dir_all, File};
use std::path::Path;
use std::str::FromStr;
use strum::{EnumProperty, IntoEnumIterator};

// ----------------------------------------------------------------------------
// Dataset Download
// ----------------------------------------------------------------------------

/// List datasets
pub async fn list_datasets(args: &cli::DatasetListArgs) -> Result<(), Report> {
    // table of name, tag, cli_version
    let mut table = utils::table::Table::new();
    table.headers = vec![
        "Name",
        "CLI Version",
        "Minimum Tag Date",
        "Maximum Tag Date",
    ]
    .into_iter()
    .map(String::from)
    .collect_vec();

    for name in Name::iter() {
        // Check if this was not the name requested by CLI args
        if let Some(args_name) = &args.name {
            if &name != args_name {
                continue;
            }
        }

        // check if this enum should be listed
        if name.get_str("list").unwrap_or("false") != "true" {
            continue;
        }

        // Extract compatibility attributes
        let compatibility = name.compatibility()?;

        let cli_version = compatibility.cli.version.unwrap_or(String::new());
        let min_date = if let Some(min_date) = compatibility.dataset.min_date {
            min_date
                .to_rfc3339_opts(SecondsFormat::Secs, true)
                .to_string()
        } else {
            String::new()
        };
        let max_date = if let Some(max_date) = compatibility.dataset.max_date {
            max_date
                .to_rfc3339_opts(SecondsFormat::Secs, true)
                .to_string()
        } else {
            "latest".to_string()
        };

        // Add to row
        let row = vec![
            name.to_string(),
            cli_version.to_string(),
            min_date.to_string(),
            max_date.to_string(),
        ];
        table.rows.push(row);
    }

    println!("\n{}", table.to_markdown()?);

    Ok(())
}

// ----------------------------------------------------------------------------
// Dataset Download
// ----------------------------------------------------------------------------

/// Download dataset
pub async fn download_dataset(args: &mut cli::DatasetDownloadArgs) -> Result<(), Report> {
    // Create the output directory if it doesn't exist
    if !args.output_dir.exists() {
        info!("Creating output directory: {:?}", &args.output_dir);
        create_dir_all(&args.output_dir)?;
    }

    // --------------------------------------------------------------------
    // Compatibility

    check_compatibility(&args.name, &args.tag)?;

    // --------------------------------------------------------------------
    // Optional Input Summary Snapshot

    let mut summary: Summary = if let Some(summary_path) = &args.summary {
        // Read in the summary file
        let reader = File::open(summary_path)
            .wrap_err_with(|| eyre!("Failed to open input summary: {summary_path:?}"))?;
        let summary: Summary = serde_json::from_reader(&reader)
            .wrap_err_with(|| eyre!("Failed to parse input summary: {summary_path:?}"))?;

        // Warn if summary conflicts with any CLI args
        if summary.name != args.name || summary.tag != args.tag {
            warn!(
                "Dataset has been changed by summary to: {} {}",
                &summary.name, &summary.tag
            );
        }
        summary
    } else {
        let mut summary = Summary::new();
        summary.name = args.name;
        summary.tag = args.tag.clone();
        summary
    };

    // --------------------------------------------------------------------
    // Download Reference

    info!("Downloading {} {} reference fasta.", &args.name, &args.tag);

    // Option #1: From Summary Snapshot
    if args.summary.is_some() {
        let ext = utils::path_to_ext(Path::new(&summary.reference.url))?;
        let decompress = ext == "zst";
        summary.reference.local_path = args.output_dir.join("reference.fasta");
        utils::download_file(
            &summary.reference.url,
            &summary.reference.local_path,
            decompress,
        )
        .await?;
    }
    // Option #2: From Dataset Tag
    else {
        summary.reference = match args.name {
            Name::SarsCov2 => sarscov2::download_reference(args).await?,
            _ => {
                return Err(eyre!(
                    "Reference download for {name} is not implemented.",
                    name = &args.name
                ))
            }
        };
    }

    // --------------------------------------------------------------------
    // Download Populations

    info!(
        "Downloading {} {} populations fasta.",
        &args.name, &args.tag
    );

    // Option #1: From Summary Snapshot
    if args.summary.is_some() {
        let ext = utils::path_to_ext(Path::new(&summary.populations.url))?;
        let decompress = ext == "zst";
        summary.populations.local_path = args.output_dir.join("populations.fasta");
        utils::download_file(
            &summary.populations.url,
            &summary.populations.local_path,
            decompress,
        )
        .await?;
    }
    // Option #2: From Dataset Tag
    else {
        summary.populations = match args.name {
            Name::SarsCov2 => sarscov2::download_populations(args).await?,
            _ => {
                return Err(eyre!(
                    "Populations download for {} is not implemented.",
                    &args.name
                ))
            }
        };
    }

    // --------------------------------------------------------------------
    // Create Annotations

    info!("Creating {} {} annotations.", &args.name, &args.tag);

    let annotations = match args.name {
        Name::SarsCov2 => sarscov2::create_annotations()?,
        _ => {
            return Err(eyre!(
                "Annotations for {} dataset is not implemented.",
                &args.name
            ))
        }
    };

    let annotations_path = args.output_dir.join("annotations.tsv");
    annotations.write(&annotations_path)?;

    // --------------------------------------------------------------------
    // Download miscellaneous files (dataset-specific)

    #[allow(clippy::single_match)]
    match args.name {
        Name::SarsCov2 => {
            // ----------------------------------------------------------------
            // Lineage Notes

            info!("Downloading {} {} lineage notes.", &args.name, &args.tag);
            // Option #1: From Summary Snapshot
            if summary.misc.contains_key("lineage_notes") {
                // Check for zst decompression
                let ext =
                    utils::path_to_ext(Path::new(&summary.misc["lineage_notes"].url))?;
                let decompress = ext == "zst";

                // Update the local path in the summary
                let mut lineage_notes_remote = summary.misc["lineage_notes"].clone();
                lineage_notes_remote.local_path =
                    args.output_dir.join("lineage_notes.txt");
                *summary.misc.get_mut("lineage_notes").unwrap() = lineage_notes_remote;

                utils::download_file(
                    &summary.misc["lineage_notes"].url,
                    &summary.misc["lineage_notes"].local_path,
                    decompress,
                )
                .await?;
            }
            // Option #2: From Dataset Tag
            else {
                let lineage_notes_remote = sarscov2::download_lineage_notes(args).await?;
                summary
                    .misc
                    .insert("lineage_notes".to_string(), lineage_notes_remote);
            }

            // ----------------------------------------------------------------
            // Alias Key

            info!("Downloading {} {} alias key.", &args.name, &args.tag);

            // Option #1: From Summary Snapshot
            if summary.misc.contains_key("alias_key") {
                let ext = utils::path_to_ext(Path::new(&summary.misc["alias_key"].url))?;
                let decompress = ext == "zst";

                // Update the local path in the summary
                let mut alias_key_remote = summary.misc["alias_key"].clone();
                alias_key_remote.local_path = args.output_dir.join("alias_key.json");
                *summary.misc.get_mut("alias_key").unwrap() = alias_key_remote;

                utils::download_file(
                    &summary.misc["alias_key"].url,
                    &summary.misc["alias_key"].local_path,
                    decompress,
                )
                .await?;
            }
            // Option #2: From Dataset Tag
            else {
                let alias_key_remote = sarscov2::download_alias_key(args).await?;
                summary
                    .misc
                    .insert("alias_key".to_string(), alias_key_remote);
            }
        }
        _ => (),
    };

    // --------------------------------------------------------------------
    // Create Phylogeny

    info!("Creating {} {} phylogeny.", &args.name, &args.tag);

    let mut phylogeny = Phylogeny::new();
    phylogeny.build_graph(args).await?;
    phylogeny.export(&args.output_dir, PhylogenyExportFormat::Dot)?;
    phylogeny.export(&args.output_dir, PhylogenyExportFormat::Json)?;

    // --------------------------------------------------------------------
    // Identify Diagnostic Mutations
    //
    // This is painfully slow, need to rethink!

    if args.diagnostic {
        let diagnostic_path = args.output_dir.join("diagnostic_mutations.tsv");
        info!(
            "Identifying {} {} diagnostic mutations.",
            &args.name, &args.tag
        );

        let mask = 0;
        let (_populations, mutations) = parse_populations(
            &summary.populations.local_path,
            &summary.reference.local_path,
            mask,
        )?;
        let diagnostic_table = get_diagnostic_mutations(&mutations, &phylogeny)?;
        diagnostic_table.write(&diagnostic_path)?;
    }

    // --------------------------------------------------------------------
    // Create Edge Cases

    info!("Creating {} {} edge cases.", &args.name, &args.tag);

    let edge_cases = match args.name {
        Name::SarsCov2 => sarscov2::create_edge_cases()?,
        _ => {
            return Err(eyre!(
                "Edge cases for {} dataset is not implemented.",
                &args.name
            ))
        }
    };
    edge_cases::export(&edge_cases, &args.output_dir, EdgeCaseExportFormat::Json)?;

    // --------------------------------------------------------------------
    // Export Summary

    info!("Exporting {} {} summary.", &args.name, &args.tag);
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

    dataset.tag = Tag::Unknown;
    dataset.name = Name::Unknown;

    let summary_path = args.dataset_dir.join("summary.json");
    if summary_path.exists() {
        info!("Loading summary: {:?}", summary_path);
        let reader = File::open(summary_path)?;
        let summary: Summary = serde_json::from_reader(&reader)?;
        dataset.name = summary.name;
        dataset.tag = summary.tag;
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
