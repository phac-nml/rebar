pub mod cli;
pub mod dataset;
pub mod export;
pub mod phylogeny;
pub mod plot;
pub mod recombination;
pub mod sequence;
pub mod utils;

use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use itertools::Itertools;
use log::{debug, info, warn};
use serde::Serialize;
use std::fs::create_dir_all;

/// Run rebar on input alignment or dataset population
pub fn run(args: cli::RunArgs) -> Result<(), Report> {
    // create output directory if it doesn't exist
    if !args.output_dir.exists() {
        info!("Creating output directory: {:?}", args.output_dir);
        create_dir_all(&args.output_dir)?;
    }

    // Collect files in dataset_dir into a dataset object
    // This mainly includes parent populations sequences
    //   and optionally a phylogenetic representation.
    let dataset = dataset::load(&args.dataset_dir, args.mask)?;

    // init a containter to hold query sequences, either dataset
    // populations or sequences from an input alignment
    let mut sequences = Vec::new();

    // ------------------------------------------------------------------------
    // Parse Input Populations

    let populations = &args.input.populations;
    if let Some(populations) = populations {
        info!("Loading query populations: {populations}");
        // split population into csv (,) parts
        let populations = populations.split(',').collect::<Vec<_>>();
        // intermediate container to hold additional populations (like descendants)
        let mut search_populations = Vec::new();

        // if population ends in '*', use phylogenetically aware mode
        // to search for it and all descendants
        for population in populations {
            if population.ends_with('*') {
                // remove last char (wildcard) with pop
                let parent = population[0..population.len() - 1].to_string();
                // get descendants of all parents
                let descendants = dataset.phylogeny.get_descendants(&parent)?;
                // add parent population and descendants
                search_populations.push(parent);
                search_populations.extend(descendants);
            } else {
                search_populations.push(population.to_string());
            }
        }

        for population in search_populations {
            // check if population is in dataset
            if dataset.populations.contains_key(&population) {
                let sequence = dataset.populations[&population].clone();
                sequences.push(sequence);
                debug!("Added population {population} to query sequences.");
            } else {
                warn!("Population {population} was not found in the dataset.");
                continue;
            }
        }
    }

    // ------------------------------------------------------------------------
    // Parse Input alignment

    let alignment = &args.input.alignment;
    if let Some(alignment) = alignment {
        info!("Loading query alignment: {:?}", alignment);
        let alignment_reader =
            fasta::Reader::from_file(alignment).expect("Unable to read alignment");

        for result in alignment_reader.records() {
            let record = result.wrap_err(format!(
                "Unable to parse alignment: {:?}",
                alignment.to_str().unwrap()
            ))?;
            let sequence = sequence::Sequence::from_record(
                record,
                Some(&dataset.reference),
                args.mask,
            )?;
            sequences.push(sequence);
        }
    }

    // ------------------------------------------------------------------------
    // Recombination Search

    info!("Running recombination search.");
    let mut best_matches = Vec::new();
    let mut recombinations = Vec::new();

    for sequence in &sequences {
        // find the best match in the dataset to the full sequence
        let best_match = dataset.search(sequence, None, None)?;

        // organize parameters for find_parents function
        let parent_search_args = recombination::parent_search::Args::new(
            &dataset,
            sequence,
            &best_match,
            &args,
        );
        // search for all recombination parents
        let (_parents, recombination) =
            recombination::parent_search::search_all(parent_search_args)?;

        best_matches.push(best_match);
        recombinations.push(recombination);
    }

    // ------------------------------------------------------------------------
    // Export Linelist (single)

    let outpath_linelist = args.output_dir.join("linelist.tsv");
    let linelist = export::LineList::create(&recombinations, &best_matches, &dataset)?;
    let linelist_table = linelist.to_table()?;
    utils::write_table(&linelist_table, &outpath_linelist, Some('\t'))?;

    // ------------------------------------------------------------------------
    // Export Barcodes (multiple, collected by recombinant)

    let outdir_barcodes = args.output_dir.join("barcodes");
    create_dir_all(&outdir_barcodes)?;

    // get unique keys of recombinants identified
    let unique_keys = recombinations
        .iter()
        .map(|rec| &rec.unique_key)
        .unique()
        .collect_vec();

    for unique_key in unique_keys {
        // filter recombinations down to just this recombinant unique_key
        let unique_rec = recombinations
            .iter()
            .filter(|rec| rec.unique_key == *unique_key)
            .cloned()
            .collect_vec();
        // combine all the sample barcode tables
        let barcode_table = recombination::combine_tables(&unique_rec)?;
        let output_barcode_table = outdir_barcodes.join(format!("{unique_key}.tsv"));
        utils::write_table(&barcode_table, &output_barcode_table, Some('\t'))?;
    }

    info!("Done.");
    Ok(())
}

/// Plot rebar output
pub fn plot(args: cli::PlotArgs) -> Result<(), Report> {
    // parse input/output paths
    let barcodes_file = &args.barcodes.barcodes_file;
    let barcodes_dir = &args.barcodes.barcodes_dir;
    let output_dir = &args.output_dir;

    // create output directory if it doesn't exist
    if !args.output_dir.exists() {
        info!("Creating output directory: {:?}", args.output_dir);
        create_dir_all(&args.output_dir)?;
    }

    // combine files from single input and directory
    let mut barcodes_files: Vec<std::path::PathBuf> = Vec::new();

    // ------------------------------------------------------------------------
    // Input File Specified

    if let Some(barcodes_file) = barcodes_file {
        barcodes_files.push(barcodes_file.clone());
    }

    // ------------------------------------------------------------------------
    // Input Directory Specified

    if let Some(barcodes_dir) = barcodes_dir {
        let files = std::fs::read_dir(barcodes_dir)?;

        for result in files {
            let file_path = result?.path();
            let file_ext = file_path.extension().ok_or_else(|| {
                eyre!("Failed to parse file extension from {file_path:?}")
            })?;
            if file_ext == "tsv" {
                barcodes_files.push(file_path.clone());
            }
        }
    }

    for barcodes_file in barcodes_files {
        info!("Plotting barcodes file: {:?}", barcodes_file);
        let output_prefix = barcodes_file.file_stem().unwrap().to_str().unwrap();
        let output_path = output_dir.join(format!("{}.png", output_prefix));
        plot::create(&barcodes_file, &args.linelist, &output_path)?;
    }

    info!("Done.");
    Ok(())
}

// ----------------------------------------------------------------------------
// Traits
// ----------------------------------------------------------------------------

pub trait ToYaml {
    fn to_yaml(&self) -> String
    where
        Self: Serialize,
    {
        serde_yaml::to_string(&self).unwrap()
    }
}
