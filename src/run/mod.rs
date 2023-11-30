use crate::cli;
use crate::dataset;
use crate::export;
use crate::recombination;

use crate::dataset::{attributes::Name, SearchResult};
use crate::recombination::Recombination;
use crate::sequence::Sequence;
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use color_eyre::Help;
use indicatif::{style::ProgressStyle, ProgressBar};
use itertools::Itertools;
use log::{debug, info, warn};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::fs::{create_dir_all, File};
use std::io::Write;

/// Run rebar on input alignment and/or dataset population(s)
pub fn run(args: &mut cli::run::Args) -> Result<(), Report> {
    // copy args for export/seralizing
    let args_export = args.clone();

    // create output directory if it doesn't exist
    if !args.output_dir.exists() {
        info!("Creating output directory: {:?}", args.output_dir);
        create_dir_all(&args.output_dir)?;
    }
    // make sure output directory is empty!
    else {
        let output_dir_is_empty = args.output_dir.read_dir()?.next().is_none();
        if !output_dir_is_empty {
            return Err(eyre!(
                "--output-dir {:?} already exists and is not empty!",
                args.output_dir
            )
            .suggestion("Please change your --output-dir to a new or empty directory."));
        }
    }

    // check how many threads are available on the system
    let default_thread_pool =
        rayon::ThreadPoolBuilder::new().build().expect("Failed to build thread pool.");
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
    let result = rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global();

    if result.is_err() {
        warn!("Failed to build global thread pool.");
    }

    // configure progress bar style
    let progress_bar_style = ProgressStyle::with_template(
        "{bar:40} {pos}/{len} ({percent}%) | Sequences / Second: {per_sec} | Elapsed: {elapsed_precise} | ETA: {eta_precise}"
    ).expect("Failed to create progress bar from template.");

    // Collect files in dataset_dir into a dataset object
    // This mainly includes parent populations sequences
    //   and optionally a phylogenetic representation.
    let mut dataset = dataset::load::dataset(&args.dataset_dir, &args.mask)?;

    // init a container to hold query sequences, dataset
    // populations and/or sequences from an input alignment
    let mut sequences = Vec::new();
    // keep track of ids we've seen to remove duplicates later
    let mut ids_seen = Vec::new();

    // ------------------------------------------------------------------------
    // Parse Input Populations
    // ------------------------------------------------------------------------

    // this step is pretty fast, don't really need a progress bar here

    if let Some(populations) = &args.input.populations {
        info!("Parsing input populations: {populations:?}");

        dataset.expand_populations(populations)?.into_iter().for_each(|p| {
            if !dataset.populations.contains_key(&p) {
                warn!("Population {p} is not in the dataset populations fasta.");
            } else {
                debug!("Adding population {p} to query sequences.");
                let mut sequence = dataset.populations.get(&p).unwrap().clone();
                sequence.id = format!("population_{}", sequence.id);
                ids_seen.push(sequence.id.clone());
                sequences.push(sequence.clone());
            }
        });
    }

    // ------------------------------------------------------------------------
    // Parse Input Alignment
    // ------------------------------------------------------------------------

    if let Some(alignment) = &args.input.alignment {
        info!("Loading query alignment: {:?}", alignment);
        let alignment_reader = fasta::Reader::from_file(alignment)
            .map_err(|e| eyre!(e))
            .wrap_err("Failed to read file: {alignment:?}")?;

        for result in alignment_reader.records() {
            let record = result.wrap_err("Unable to parse alignment: {alignment:?}")?;
            let sequence =
                Sequence::from_record(record, Some(&dataset.reference), &args.mask)?;

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
    // Parse and expand input parents

    if let Some(parents) = &args.parents {
        info!("Parsing input parents: {:?}", &parents);
        args.parents = Some(dataset.expand_populations(parents)?);
    }

    // ------------------------------------------------------------------------
    // Dataset Knockout
    // ------------------------------------------------------------------------

    if let Some(knockout) = &args.knockout {
        info!("Performing dataset knockout: {knockout:?}");

        let mut expanded_knockout = Vec::new();

        for p in knockout {
            // Replace the wildcard, knockout will always include descendants
            let p = p.replace('*', "");

            // Identify descendants
            let exclude_populations = if dataset.phylogeny.is_empty() {
                vec![p.to_string()]
            } else {
                dataset.phylogeny.get_descendants(&p)?
            };

            // remove from populations
            debug!("Removing {p}* from the populations fasta.");
            dataset.populations.retain(|id, _| !exclude_populations.contains(id));

            // remove from mutations
            debug!("Removing {p}* from the mutations.");
            dataset.mutations = dataset
                .mutations
                .into_iter()
                .filter_map(|(sub, mut populations)| {
                    populations.retain(|p| !exclude_populations.contains(p));
                    if populations.is_empty() {
                        None
                    } else {
                        Some((sub, populations))
                    }
                })
                .collect();

            // remove from phylogeny
            if !dataset.phylogeny.is_empty() {
                debug!("Removing {p}* from the phylogeny.");
                dataset.phylogeny = dataset.phylogeny.prune(&p)?;
            }

            expanded_knockout.extend(exclude_populations);
        }

        args.knockout = Some(expanded_knockout);
    }

    // ------------------------------------------------------------------------
    // Recombination Search
    // ------------------------------------------------------------------------

    info!("Running recombination search.");

    // this step is the slowest, use progress bar and parallel threads
    let progress_bar = ProgressBar::new(sequences.len() as u64);
    progress_bar.set_style(progress_bar_style);

    // adjust search populations based on args.parents and args.knockout
    let mut parent_search_populations = dataset.populations.keys().collect_vec();
    // if args.parents supplied on the CLI
    if let Some(populations) = &args.parents {
        parent_search_populations.retain(|pop| populations.contains(pop))
    }
    // if args.knockout supplied on the CLI
    if let Some(populations) = &args.knockout {
        parent_search_populations.retain(|pop| !populations.contains(pop))
    }

    // Search for the best match and recombination parents for each sequence.
    // This loop/closure is structured weirdly for rayon compatability, and the
    // fact that we need to return multiple types of objects
    let results: Vec<(SearchResult, Recombination)> = sequences
        .par_iter()
        .map(|sequence| {
            // initialize with default results, regardless of whether our
            // searches "succeed", we're going to return standardized data
            // structures to build our exports upon (ex. linelist columns)
            // which will include the "negative" results
            let mut best_match = SearchResult::new(sequence);
            let mut recombination = Recombination::new(sequence);

            // ------------------------------------------------------------------------
            // Best Match (Consensus)
            //
            // search for the best match in the dataset to this sequence.
            // this will represent the consensus population call.

            debug!("Identifying best match (consensus population).");
            let search_result = dataset.search(sequence, None, None);

            // if we found a match, proceed with recombinant search
            if let Ok(search_result) = search_result {
                best_match = search_result;

                debug!("Searching for recombination parents.");
                let parent_search = recombination::search::all_parents(
                    sequence,
                    &dataset,
                    &mut best_match,
                    &parent_search_populations,
                    args,
                );
                match parent_search {
                    Ok(search_result) => recombination = search_result,
                    Err(e) => debug!("Parent search did not succeed. {e}"),
                }
            }
            // what to do if not a single population matched?
            else {
                // temporary handling for root population B
                if dataset.name == Name::SarsCov2 {
                    if sequence.id == "population_B" {
                        best_match.consensus_population = "B".to_string();
                    }
                } else {
                    debug!("No matches found.");
                }
            }

            progress_bar.inc(1);

            (best_match, recombination)
        })
        .collect();

    progress_bar.finish();

    // ------------------------------------------------------------------------
    // Export CLI args

    let outpath_args = args.output_dir.join("run_args.json");
    info!("Exporting CLI Run Args: {outpath_args:?}");
    // create output file
    let mut file = File::create(&outpath_args)
        .wrap_err_with(|| format!("Failed to create file: {outpath_args:?}"))?;

    // parse to string
    let output = serde_json::to_string_pretty(&args_export)
        .wrap_err_with(|| "Failed to parse mutations.".to_string())?;

    // write to file
    file.write_all(format!("{}\n", output).as_bytes())
        .wrap_err_with(|| format!("Failed to write file: {outpath_args:?}"))?;

    // ------------------------------------------------------------------------
    // Export Linelist (single)

    let outpath_linelist = args.output_dir.join("linelist.tsv");
    info!("Exporting linelist: {outpath_linelist:?}");

    let linelist_table = export::linelist(&results, &dataset)?;
    //let linelist_table = export::linelist(&best_matches, &recombinations, &dataset)?;
    linelist_table.write(&outpath_linelist)?;

    // ------------------------------------------------------------------------
    // Export Barcodes (multiple, collected by recombinant)

    let outdir_barcodes = args.output_dir.join("barcodes");

    create_dir_all(&outdir_barcodes)?;

    // get unique keys of recombinants identified
    let unique_keys = results
        .iter()
        .filter_map(|(_b, r)| (*r.unique_key != String::new()).then_some(&r.unique_key))
        .unique()
        .collect_vec();

    if unique_keys.is_empty() {
        warn!("No recombination detected, no barcodes will be outputted.");
    } else {
        info!("Exporting recombination barcodes: {outdir_barcodes:?}");
    }

    for unique_key in unique_keys {
        // filter recombinations down to just this recombinant unique_key
        let unique_rec = results
            .iter()
            .filter_map(|(_b, r)| (r.unique_key == *unique_key).then_some(r))
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
