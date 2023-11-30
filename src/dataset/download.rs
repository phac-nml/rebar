use crate::cli;
use crate::dataset;
use crate::dataset::attributes::{check_compatibility, Name, Summary};
use crate::dataset::{sarscov2, toy1};
use crate::{utils, utils::remote_file::RemoteFile};
use color_eyre::eyre::{Report, Result};
use itertools::Itertools;
use log::{info, warn};
use std::fs::create_dir_all;
use std::path::Path;

/// Download dataset
pub async fn dataset(args: &mut cli::dataset::download::Args) -> Result<(), Report> {
    info!("Downloading dataset: {} {}", &args.name, &args.tag);

    // --------------------------------------------------------------------
    // Optional Input Summary Snapshot

    let mut summary: Summary = if let Some(summary_path) = &args.summary {
        info!("Importing summary: {summary_path:?}");
        let summary = Summary::read(summary_path)?;

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
    // Compatibility Check

    check_compatibility(&args.name, &args.tag)?;

    // Warn if the directory already exists
    if !args.output_dir.exists() {
        info!("Creating output directory: {:?}", &args.output_dir);
        create_dir_all(&args.output_dir)?;
    } else {
        warn!(
            "--output-dir {:?} already exists, proceed with caution!",
            args.output_dir
        );
    }

    // --------------------------------------------------------------------
    // Reference

    let output_path = args.output_dir.join("reference.fasta");
    info!("Downloading reference: {output_path:?}");

    summary.reference = if args.summary.is_some() {
        snapshot(&summary.reference, &output_path).await?
    } else {
        match args.name {
            Name::SarsCov2 => {
                sarscov2::download::reference(&args.tag, &output_path).await?
            }
            Name::Toy1 => toy1::download::reference(&args.tag, &output_path)?,
            _ => todo!(),
        }
    };

    // --------------------------------------------------------------------
    // Populations

    let output_path = args.output_dir.join("populations.fasta");
    info!("Downloading populations: {output_path:?}");

    summary.populations = if args.summary.is_some() {
        dataset::download::snapshot(&summary.populations, &output_path).await?
    } else {
        match args.name {
            Name::SarsCov2 => {
                sarscov2::download::populations(&args.tag, &output_path).await?
            }
            Name::Toy1 => toy1::download::populations(&args.tag, &output_path)?,
            _ => todo!(),
        }
    };

    // --------------------------------------------------------------------
    // Annotations

    let output_path = args.output_dir.join("annotations.tsv");
    info!("Creating annotations: {output_path:?}");

    let annotations = match args.name {
        Name::SarsCov2 => sarscov2::annotations::build()?,
        Name::Toy1 => toy1::annotations::build()?,
        _ => todo!(),
    };
    annotations.write(&output_path)?;

    // --------------------------------------------------------------------
    // Graph (Phylogeny)

    let output_path = args.output_dir.join("phylogeny.json");
    info!("Building phylogeny: {output_path:?}");

    let phylogeny = match args.name {
        Name::SarsCov2 => {
            sarscov2::phylogeny::build(&mut summary, &args.output_dir).await?
        }
        Name::Toy1 => toy1::phylogeny::build()?,
        _ => todo!(),
    };
    phylogeny.write(&output_path)?;
    // Also write as .dot file for graphviz visualization.
    let output_path = args.output_dir.join("phylogeny.dot");
    info!("Exporting graphviz phylogeny: {output_path:?}");
    phylogeny.write(&output_path)?;

    // --------------------------------------------------------------------
    // Export Mutations

    let output_path = args.output_dir.join("mutations.json");
    info!("Mapping mutations to populations: {output_path:?}");
    let mask = vec![0, 0];
    let (_populations, mutations) = dataset::load::parse_populations(
        &summary.populations.local_path,
        &summary.reference.local_path,
        &mask,
    )?;
    dataset::write_mutations(&mutations, &output_path)?;

    // --------------------------------------------------------------------
    // Create Edge Cases
    //
    // Edge cases are simply a vector of the CLI Run Args (cli::run::Args)
    // customized to particular recombinants.

    let output_path = args.output_dir.join("edge_cases.json");
    info!("Creating edge cases: {output_path:?}");

    let mut edge_cases = match args.name {
        Name::SarsCov2 => dataset::sarscov2::edge_cases::default()?,
        Name::Toy1 => dataset::toy1::edge_cases::default()?,
        _ => todo!(),
    };
    let manual_populations =
        edge_cases.iter().filter_map(|e| e.population.clone()).collect_vec();

    let problematic_recombinants = phylogeny.get_problematic_recombinants()?;
    for recombinant in problematic_recombinants {
        let parents = phylogeny.get_parents(&recombinant)?;
        warn!("Recombinant {recombinant} is problematic. Parents are not sister taxa: {parents:?}");
        if manual_populations.contains(&recombinant) {
            warn!("Manual edge case exists: {recombinant:?}");
        } else {
            warn!("Creating auto edge case: {recombinant:?}");
            let edge_case = cli::run::Args {
                population: Some(recombinant),
                parents: Some(parents),
                ..Default::default()
            };
            edge_cases.push(edge_case);
        }
    }

    // Reminder, we use the module write  method, not the struct method,
    // because this is a vector of arguments we need to serialize.
    cli::run::Args::write(&edge_cases, &output_path)?;

    // --------------------------------------------------------------------
    // Export Summary

    let output_path = args.output_dir.join("summary.json");
    info!("Exporting summary: {output_path:?}");
    summary.write(&output_path)?;

    // --------------------------------------------------------------------
    // Finish

    info!("Done.");
    Ok(())
}

/// Download remote file from a summary snapshot.
pub async fn snapshot(
    snapshot: &RemoteFile,
    output_path: &Path,
) -> Result<RemoteFile, Report> {
    // Check extension for decompression
    let ext = utils::path_to_ext(Path::new(&snapshot.url))?;
    let decompress = ext == "zst";

    // Update the local path to the desired output
    let mut remote_file = snapshot.clone();
    remote_file.local_path = output_path.to_path_buf();

    // Download the file
    utils::download_file(&snapshot.url, output_path, decompress).await?;

    Ok(remote_file)
}
