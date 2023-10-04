pub mod cli;
pub mod dataset;
pub mod export;
pub mod phylogeny;
pub mod plot;
pub mod recombination;
pub mod sequence;
pub mod utils;

use crate::dataset::attributes::Name;
use crate::dataset::SearchResult;
use crate::recombination::{search, Recombination};
use crate::sequence::Sequence;
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use indicatif::{style::ProgressStyle, ProgressBar};
use itertools::Itertools;
use log::{debug, info, warn};
use rayon::prelude::*;
use std::fs::create_dir_all;

/// Download rebar dataset
pub async fn download_dataset(
    args: &mut cli::dataset::download::Args,
) -> Result<(), Report> {
    dataset::download::dataset(args).await?;
    Ok(())
}

/// List rebar datasets
pub async fn list_datasets(args: &cli::dataset::list::Args) -> Result<(), Report> {
    dataset::list::datasets(args).await?;
    Ok(())
}

/// Run rebar on input alignment and/or dataset population
pub fn run(args: cli::run::Args) -> Result<(), Report> {
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
    let mut dataset = dataset::load::dataset(&args)?;

    // init a container to hold query sequences, dataset
    // populations and/or sequences from an input alignment
    let mut sequences = Vec::new();
    // keep track of ids we've seen to remove duplicates later
    let mut ids_seen = Vec::new();

    // ------------------------------------------------------------------------
    // Parse Input Populations

    // this step is pretty fast, don't really need a progress bar here

    if let Some(populations) = &args.input.populations {
        info!("Loading query populations: {populations:?}");

        // Parse populations, and expand '*' to get descendants
        let search_populations = populations
            .into_iter()
            .map(|p| {
                // if population is '*', use all populations in dataset
                if p == "*" {
                    Ok(dataset.populations.keys().cloned().collect_vec())
                }
                // if population ends with '*' expand descendants
                else if p.ends_with('*') {
                    let p = p.replace("*", "");
                    dataset.phylogeny.get_descendants(&p)
                } 
                // simple population name, that is in the dataset
                else if dataset.populations.contains_key(p) {
                    Ok(vec![p.to_string()])
                }
                else {
                    Err(eyre!("{p} is not present in the dataset."))
                }
            })
            // flatten and handle the `Result` layer
            .collect::<Result<Vec<_>, Report>>()?
            .into_iter()
            // flatten and handle the `Vec` layer
            .flatten()
            .unique()
            .collect_vec();

        (ids_seen, sequences) = search_populations
            .into_iter()
            .map(|p| {
                debug!("Adding population {p} to query sequences.");
                let mut sequence = dataset.populations[&p].clone();
                sequence.id = format!("population_{}", sequence.id).to_string();
                (sequence.id.clone(), sequence)
            })
            .unzip();
    }

    // ------------------------------------------------------------------------
    // Parse Input alignment

    if let Some(alignment) = &args.input.alignment {
        info!("Loading query alignment: {:?}", alignment);
        let alignment_reader = fasta::Reader::from_file(alignment)
            .map_err(|e| eyre!(e))
            .wrap_err("Failed to read file: {alignment:?}")?;

        for result in alignment_reader.records() {
            let record = result.wrap_err("Unable to parse alignment: {alignment:?}")?;
            let sequence = Sequence::from_record(record, Some(&dataset.reference), &args.mask)?;

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
    // Dataset Knockout

    if let Some(knockout) = &args.knockout {
        info!("Performing dataset knockout: {knockout:?}");

        //let mut knockout_populations = vec![];
    
        for p in knockout {
            let p = p.replace("*", "");

            // remove from populations
            debug!("Removing {p} from the populations fasta.");
            let exclude_populations = dataset.phylogeny.get_descendants(&p)?;
            dataset.populations.retain(|id, _| !exclude_populations.contains(id));

            // remove from phylogeny
            debug!("Removing {p} from the phylogeny.");
            dataset.phylogeny = dataset.phylogeny.prune(&p)?;

            // remove from mutations
            debug!("Removing {p} from the mutations.");
            dataset.mutations = dataset
                .mutations
                .into_iter()
                .filter_map(|(sub, mut populations)| {
                    populations.retain(|p| !exclude_populations.contains(&p));
                    if populations.is_empty() { None } 
                    else { Some((sub, populations))}
                })
                .collect();
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
            let mut best_match = SearchResult::new(sequence);
            let mut recombination = Recombination::new(sequence);

            // search for the best match in the dataset to this sequence
            let search_result = dataset.search(sequence, None, None);

            if let Ok(search_result) = search_result {
                // use the successful search as the best_match
                best_match = search_result;
                // search for all recombination parents (primary and then secondary)
                let parent_search =
                    search::all_parents(&sequence, &dataset, &mut best_match, &args);
                // if the search was successful, unpack the results or issue a warning
                match parent_search {
                    Ok(search_result) => (_, recombination) = search_result,
                    Err(e) => warn!("{e:?}")
                }
            }
            // what to do if not a single population matched?
            else {
                // temporary handling for root population B
                if dataset.name == Name::SarsCov2 {
                    if sequence.id == "population_B".to_string() {
                        best_match.consensus_population = "B".to_string();
                    }
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
        warn!("No recombination detected, no barcodes will be outputted.");
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
pub fn plot(args: cli::plot::Args) -> Result<(), Report> {
    // ------------------------------------------------------------------------
    // Parse Args

    let dataset_dir = &args.dataset_dir;
    if !dataset_dir.exists() {
        return Err(eyre!("--dataset-dir {dataset_dir:?} does not exist."));
    }
    let output_dir = &args.output_dir;
    if !output_dir.exists() {
        return Err(eyre!("--output-dir {output_dir:?} does not exist."));
    }

    // ------------------------------------------------------------------------
    // Check Mandatory Paths

    let linelist = &output_dir.join("linelist.tsv");
    if !linelist.exists() {
        return Err(eyre!(
            "Linelist file {linelist:?} does not exist in --output-dir {output_dir:?}."
        ));
    }
    let barcodes_dir = &output_dir.join("barcodes");
    if !linelist.exists() {
        return Err(eyre!("Barcodes directory {barcodes_dir:?} does not exist in --output-dir {output_dir:?}."));
    }
    let annotations = &dataset_dir.join("annotations.tsv");
    if !annotations.exists() {
        return Err(eyre!(
            "Annotations {annotations:?} do not exist in --dataset-dir {dataset_dir:?}."
        ));
    }

    // create plot directory if it doesn't exist
    let plot_dir = args.plot_dir.unwrap_or(output_dir.join("plots"));
    if !plot_dir.exists() {
        info!("Creating plot directory: {plot_dir:?}");
        create_dir_all(&plot_dir)?;
    }

    // ------------------------------------------------------------------------
    // List of Barcodes Files

    let mut barcodes_files: Vec<std::path::PathBuf> = Vec::new();

    // Input File Specified
    let barcodes_file = &args.barcodes_file;
    if let Some(barcodes_file) = barcodes_file {
        if !barcodes_file.exists() {
            return Err(eyre!("Barcodes file {barcodes_file:?} does not exist."));
        }
        barcodes_files.push(barcodes_file.clone());
    }
    // Otherwise, use all files in barcodes_dir
    else {
        let files = std::fs::read_dir(barcodes_dir)?;
        for result in files {
            let file_path = result?.path();
            let file_ext = file_path.extension().unwrap_or(std::ffi::OsStr::new(""));
            if file_ext == "tsv" {
                barcodes_files.push(file_path.clone());
            } else {
                warn!("Skipping barcodes file with unknown extension: {file_path:?}")
            }
        }
    }

    // ------------------------------------------------------------------------
    // Plot Each Barcodes

    for barcodes_file in barcodes_files {
        info!("Plotting barcodes file: {:?}", barcodes_file);
        let output_prefix = barcodes_file.file_stem().unwrap().to_str().unwrap();
        let output_path = plot_dir.join(format!("{}.png", output_prefix));
        let result =
            plot::create(&barcodes_file, linelist, Some(annotations), &output_path);
        match result {
            Ok(_) => info!("Plotting success."),
            Err(e) => warn!("Plotting failure. The following error was encountered but ignored: {e:?}"),
        }
    }

    info!("Done.");
    Ok(())
}
