// System paths
// use std::path::Path;
//use tempfile::Builder;
// Logging
//use log::info;
use clap::Parser;
use color_eyre::eyre::{Report, Result};
use log::{debug, info};
use rebar::cli::verbosity::Verbosity;
use rebar::cli::{Cli, Command};
use rebar::dataset::Dataset;
use rebar::query::Query;
use rebar::traits::ToYaml;
use std::env;
use std::str::FromStr;

fn setup() -> Result<(), Report> {
    color_eyre::install()?;

    let verbosity = Verbosity::from_str("debug").unwrap();

    // Set default logging level if RUST_LOG is not set.
    //let log_verbosity = "debug";
    if env::var("RUST_LOG").is_err() {
        env::set_var("RUST_LOG", verbosity.to_string())
    }

    env_logger::init();

    Ok(())
}

#[tokio::main]
async fn main() -> Result<()> {
    // Misc setup actions like logging
    setup().unwrap();

    // Parse CLI parameters
    let args = Cli::parse();

    match args.command {
        // Download a dataset
        Command::Dataset {
            name,
            tag,
            output_dir,
            ..
        } => {
            Dataset::download(&name, &tag, &output_dir).await?;
        }
        // Run on input alignment
        Command::Run {
            alignment,
            dataset_dir,
            mask,
            ..
        } => {
            // Load dataset
            let dataset = Dataset::load(&dataset_dir, mask)?;
            // Load the query alignment
            // let query = Query::load(alignment, &dataset, mask)?;
            // info!("Identifying consensus and parent populations.");
            // for (id, sequence) in query.sequences {

            //     if id != "XBB.1.16" { continue }
            //     debug!("sequence: {id}");

            //     // consensus population
            //     let exclude_populations = None;
            //     let best_match = dataset.find_best_match(&sequence, exclude_populations)?;
            //     debug!("\n  {}", best_match.to_yaml().replace("\n", "\n  "));

            //     // let mut exclude_populations = vec!(best_match.consensus_population);
            //     // exclude_populations.push(String::from("XBB.1.16.1"));
            //     // exclude_populations.push(String::from("XBB.1.16.3"));
            //     // let parent_match = dataset.find_best_match(&sequence, Some(exclude_populations))?;
            //     // debug!("  parent_1:");
            //     // debug!("\n    {}", parent_match.to_yaml().replace("\n", "\n    "));
            //     // parent populations
            // }
        }
    }

    // ------------------------------------------------------------------------
    // Dataset Creation

    //dataset.populations = dataset.download_sequences(&dataset_outdir)?;
    //let dataset = Dataset::create(dataset_name.to_string(), dataset_tag.to_string(), mask);

    // Sequences
    // dataset.populations.set_sequences(&reference_path, &populations_path, &mask).unwrap();
    // dataset.populations.set_mutations().unwrap();

    // // Phylogeny
    // info!("Preparing dataset phylogeny: {}", &phylogeny_path.display());
    // let phylogeny
    // dataset.phylogeny.build_graph(&dataset_name, &dataset_tag, &dataset_dir).expect("Failed to build phylogeny.");
    // dataset.phylogeny.export_graph(&dataset_dir).expect("Failed to export phylogeny.");

    // // Sequences
    // info!("Preparing dataset sequences: {}", &populations_path.display());
    // dataset.populations.set_sequences(&reference_path, &populations_path, &mask).unwrap();
    // dataset.populations.set_mutations().unwrap();

    // // ------------------------------------------------------------------------
    // // Run

    // info!("Importing query sequences: {}", &sequences_path.display());
    // let mut query = Sequences::new();
    // query.set_sequences(&reference_path, &sequences_path, &mask).unwrap();
    // query.summarise_barcodes(&dataset).unwrap();

    // 1. Sumarise the barcodes (conflicts, support, total), etc.
    // 2. Assign consensus population (Barcode Match)
    // 3. Find recombination parents (Barcode Match)

    Ok(())
}
