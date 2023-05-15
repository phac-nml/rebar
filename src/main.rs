use std::env;

// System paths
// use std::path::Path;
//use tempfile::Builder;

// Logging
//use log::info;
use color_eyre::eyre::Result;
use rebar::cli::log::LogVerbosity;

// Dataset
use rebar::dataset::Dataset;

// Traits
use std::str::FromStr;

#[tokio::main]
async fn main() -> Result<()> {
    let log_verbosity = LogVerbosity::from_str("debug").unwrap();

    // Set default logging level if RUST_LOG is not set.
    //let log_verbosity = "debug";
    if env::var("RUST_LOG").is_err() {
        env::set_var("RUST_LOG", log_verbosity.to_string())
    }

    env_logger::init();

    color_eyre::install()?;

    let dataset_name = "sars-cov-2";
    let dataset_version = "nightly";
    // let mask = 200;

    // ------------------------------------------------------------------------
    // Dataset Creation

    let dataset =
        Dataset::new(dataset_name.to_string(), dataset_version.to_string())?;
    println!("{}", dataset);
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
