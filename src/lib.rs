pub mod cli;
pub mod dataset;
pub mod export;
pub mod phylogeny;
pub mod plot;
pub mod recombination;
pub mod sequence;
pub mod utils;

use crate::dataset::SearchResult;
use crate::recombination::{search, Recombination};
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use indicatif::{style::ProgressStyle, ProgressBar};
use itertools::Itertools;
use log::{debug, info, warn};
use rayon::prelude::*;
use std::fs::create_dir_all;

/// Run rebar on input alignment and/or dataset population
pub fn run(args: cli::RunArgs) -> Result<(), Report> {
    // check how many threads are available on the system
    let default_thread_pool = rayon::ThreadPoolBuilder::new().build().unwrap();
    info!(
        "Number of threads available: {}",
        default_thread_pool.current_num_threads()
    );

    // warn the user if they requested more than their system has available
    // if so, default to the system threads
    let mut num_threads = args.threads;
    if args.threads > default_thread_pool.current_num_threads() {
        warn!(
            "Requested --threads {} is greater than the available threads.",
            args.threads
        );
        num_threads = default_thread_pool.current_num_threads();
    }

    // configure the global thread pool
    info!("Using {} thread(s).", num_threads);
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();

    // configure progress bar style
    let progress_bar_style = ProgressStyle::with_template(
        "{bar:40} {pos}/{len} ({percent}%) | Sequences / Second: {per_sec} | Elapsed: {elapsed_precise} | ETA: {eta_precise}"
    ).unwrap();

    // create output directory if it doesn't exist
    if !args.output_dir.exists() {
        info!("Creating output directory: {:?}", args.output_dir);
        create_dir_all(&args.output_dir)?;
    }

    // Collect files in dataset_dir into a dataset object
    // This mainly includes parent populations sequences
    //   and optionally a phylogenetic representation.
    let dataset = dataset::io::load_dataset(&args)?;

    // init a container to hold query sequences, dataset
    // populations and/or sequences from an input alignment
    let mut sequences = Vec::new();
    // keep track of ids we've seen to remove duplicates later
    let mut ids_seen = Vec::new();

    // ------------------------------------------------------------------------
    // Parse Input Populations

    // this step is pretty fast, don't really need a progress bar here

    let populations = &args.input.populations;
    if let Some(populations) = populations {
        info!("Loading query populations: {populations}");
        // split population into csv (,) parts
        let populations = populations.split(',').collect::<Vec<_>>();
        // intermediate container to hold additional populations (like descendants)
        let mut search_populations = Vec::new();

        for population in populations {
            // if population is '*', use all populations in dataset
            if population == "*" {
                search_populations = dataset.populations.keys().cloned().collect_vec();
                // break from loop, we're using all populations
                break;
            }
            // if population ends in '*', use phylogenetically aware mode
            // to search for it and all descendants
            if population.ends_with('*') {
                // remove last char (wildcard) with pop
                let parent = population[0..population.len() - 1].to_string();
                // get descendants of this population
                let descendants = dataset.phylogeny.get_descendants(&parent)?;
                // add population and its descendants to our search list
                search_populations.push(parent);
                search_populations.extend(descendants);
            } else {
                // no phylogenetically aware mode, just add this
                // population to our search list
                search_populations.push(population.to_string());
            }
        }

        // de-duplicate search populations
        let search_populations = search_populations.into_iter().unique().collect_vec();

        for population in search_populations {
            // check if population is in dataset
            if dataset.populations.contains_key(&population) {
                let mut sequence = dataset.populations[&population].clone();
                // add prefix 'population_' to differentiate from dataset itself
                sequence.id = format!("population_{}", sequence.id).to_string();
                ids_seen.push(sequence.id.clone());
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

            // check for duplicates
            if ids_seen.contains(&sequence.id) {
                warn!(
                    "Sequence {} is duplicated, retaining first one.",
                    sequence.id
                );
                continue;
            } else {
                ids_seen.push(sequence.id.clone());
                sequences.push(sequence);
            }
        }
    }

    // ------------------------------------------------------------------------
    // Recombination Search

    info!("Running recombination search.");

    // this step is the slowest, use progress bar and parallel threads
    let progress_bar = ProgressBar::new(sequences.len() as u64);
    progress_bar.set_style(progress_bar_style);

    // Search for the best matches and recombination parents for each sequence.
    // This loop/closure is structured weirdly for rayon compatability, and the
    // fact that we need to return multiple types of objects
    let results: Vec<(SearchResult, Recombination)> = sequences
        .par_iter()
        .map(|sequence| {
            // initialize with default results, regardless of whether our
            // searches "succeed", we're going to return standardized data
            // structures to build our exports upon (ex. linelist columns)
            let mut best_match = SearchResult::default();
            let mut recombination = Recombination::new(sequence);

            // search for the best match in the dataset to this sequence
            let search_result = dataset.search(sequence, None, None);
            if let Ok(search_result) = search_result {
                // use the successful search as the best_match
                best_match = search_result;
                // setup a search for the recombination parents
                let mut parent_search_args =
                    search::Args::new(&dataset, sequence, &best_match, &args);
                let parent_search = search::all_parents(&mut parent_search_args);
                // if the search was successful, unzip the results
                if let Ok(parent_search) = parent_search {
                    (_, recombination) = parent_search;
                }
            }

            progress_bar.inc(1);

            (best_match, recombination)
        })
        .collect();

    progress_bar.finish();

    // unpack our search results, should be able to do this with unzip, but I
    // keep getting errors related to implementation of Default and Extend

    let mut best_matches = Vec::new();
    let mut recombinations = Vec::new();

    for result in results {
        best_matches.push(result.0);
        recombinations.push(result.1);
    }

    // ------------------------------------------------------------------------
    // Export Linelist (single)

    let outpath_linelist = args.output_dir.join("linelist.tsv");
    info!("Exporting linelist: {outpath_linelist:?}");

    let linelist_table = export::linelist(&recombinations, &best_matches, &dataset)?;
    linelist_table.write(&outpath_linelist)?;

    // ------------------------------------------------------------------------
    // Export Barcodes (multiple, collected by recombinant)

    let outdir_barcodes = args.output_dir.join("barcodes");

    create_dir_all(&outdir_barcodes)?;

    // get unique keys of recombinants identified
    let unique_keys = recombinations
        .iter()
        .filter(|rec| *rec.unique_key != String::new())
        .map(|rec| &rec.unique_key)
        .unique()
        .collect_vec();

    if unique_keys.is_empty() {
        warn!("No recombination detected, no barcodes will be output");
    } else {
        info!("Exporting recombination barcodes: {outdir_barcodes:?}");
    }

    for unique_key in unique_keys {
        // filter recombinations down to just this recombinant unique_key
        let unique_rec = recombinations
            .iter()
            .filter(|rec| rec.unique_key == *unique_key)
            .cloned()
            .collect_vec();
        // combine all the sample barcode tables
        let barcode_table =
            recombination::combine_tables(&unique_rec, &dataset.reference)?;
        let barcode_table_path = outdir_barcodes.join(format!("{unique_key}.tsv"));
        barcode_table.write(&barcode_table_path)?;
    }

    info!("Done.");
    Ok(())
}

/// Plot rebar output
pub fn plot(args: cli::PlotArgs) -> Result<(), Report> {
    // parse mandatory paths
    let linelist = &args.linelist;
    let output_dir = &args.output_dir;

    // parse optional paths
    let barcodes_file = &args.barcodes.barcodes_file;
    let barcodes_dir = &args.barcodes.barcodes_dir;
    // Convert annotations from Option<PathBuf> to Option<&Path> for plot::create
    let annotations = args.annotations.as_deref();

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
            // error we're debugging where just a ".tsv" file is produced
            if file_path == barcodes_dir.join(".tsv") {
                continue;
            }
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
        plot::create(&barcodes_file, linelist, annotations, &output_path)?;
    }

    info!("Done.");
    Ok(())
}
