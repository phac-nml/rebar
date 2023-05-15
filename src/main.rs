use std::env;

// System paths
use std::path::Path;
use tempfile::Builder;

// Logging
use log::info;
use rebar::cli::log::LogVerbosity;

// Dataset
use rebar::dataset::Dataset;

// Traits
use std::str::FromStr;
use rebar::traits::ToYaml;

#[tokio::main]
async fn main() -> Result<(), std::io::Error> {
    let log_verbosity = LogVerbosity::from_str("debug").unwrap();

    // Set default logging level if RUST_LOG is not set.
    //let log_verbosity = "debug";
    if env::var("RUST_LOG").is_err() {
        env::set_var("RUST_LOG", log_verbosity.to_string())
    }

    env_logger::init();

    let dataset_name = "sars-cov-2";
    let dataset_tag = "nightly";
    let dataset_dir = Path::new("dataset/sars-cov-2/nightly");

    let mask = 200;

    let sequences_path = Path::new("data/XBB.1.16.fasta");
    let reference_path = dataset_dir.join("reference.fasta");
    let populations_path = dataset_dir.join("sequences.fasta");
    let phylogeny_path = dataset_dir.join("graph.dot");


    // let tmp_dir = Builder::new().prefix("example").tempdir()?;
    // println!("{:?}", tmp_dir);
    // let target = "https://www.rust-lang.org/logos/rust-logo-512x512.png";
    // let response = reqwest::get(target).await?;

    // ------------------------------------------------------------------------
    // Dataset


    let mut dataset = Dataset::new(dataset_name.to_string(), dataset_tag.to_string());

    // // Sequences
    // info!("Preparing dataset sequences: {}", &populations_path.display());
    // dataset.populations.set_sequences(&reference_path, &populations_path, &mask).unwrap();
    // dataset.populations.set_mutations().unwrap();
    
    println!("{}", dataset.to_yaml());
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
