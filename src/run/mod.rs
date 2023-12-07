use crate::cli;
use crate::dataset;
use crate::export;
use crate::recombination;

use crate::dataset::{attributes::Name, Dataset, SearchResult};
use crate::recombination::Recombination;
use crate::sequence::Sequence;
use crate::utils;

use color_eyre::eyre::{Report, Result, WrapErr};
use indicatif::{style::ProgressStyle, ProgressBar};
use itertools::Itertools;
use log::{debug, info, warn};
use noodles::fasta;
use rayon::iter::{ParallelBridge, ParallelIterator};
use std::fs::{create_dir_all, File, OpenOptions};
use std::io::{BufReader, Write};
use std::sync::Mutex;

/// Run rebar on input alignment and/or dataset population(s)
pub fn run(args: &mut cli::run::Args) -> Result<(), Report> {
    // Warn if the directory already exists
    if !args.output_dir.exists() {
        info!("Creating output directory: {:?}", &args.output_dir);
        create_dir_all(&args.output_dir)?;
    } else {
        warn!("Output directory already exists: {:?}", args.output_dir);
    }

    // ------------------------------------------------------------------------
    // Export CLI args
    // ------------------------------------------------------------------------

    let path = args.output_dir.join("run_args.json");
    info!("Exporting CLI Run Args: {path:?}");

    // create output file
    let mut file = File::create(&path)
        .wrap_err_with(|| format!("Failed to create file: {path:?}"))?;

    // parse to string
    let output =
        serde_json::to_string_pretty(args).wrap_err("Failed to CLI Run Args.")?;

    // write to file
    file.write_all(format!("{}\n", output).as_bytes())
        .wrap_err_with(|| format!("Failed to write file: {path:?}"))?;

    // ------------------------------------------------------------------------
    // Threads

    // check how many threads are available on the system
    let default_thread_pool = rayon::ThreadPoolBuilder::new()
        .build()
        .wrap_err("Failed to build thread pool.")?;
    let available_threads = default_thread_pool.current_num_threads();
    info!("Number of threads available: {available_threads}");

    // warn the user if they requested more than their system has available
    // if so, default to the system threads
    let mut num_threads = args.threads;
    if num_threads > available_threads {
        warn!("--threads {num_threads} is greater than the available threads.");
        num_threads = available_threads;
    }

    // configure the global thread pool
    info!("Using {} thread(s).", num_threads);
    let result = rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global();

    // the global thread pool might fail if it has already been initialized,
    // we've seen this integration tests and unittests
    if result.is_err() {
        warn!("Failed to build global thread pool.");
    }

    // ------------------------------------------------------------------------
    // Dataset
    // ------------------------------------------------------------------------

    // Collect files in dataset_dir into a dataset object
    // This mainly includes parent populations sequences
    //   and optionally a phylogensetic representation.
    let mut dataset = dataset::load::dataset(&args.dataset_dir, &args.mask)?;

    // adjust search populations based on args.parents and args.knockout
    // by default, use all populations in dataset
    let mut parent_search_populations = dataset.populations.keys().cloned().collect_vec();

    // ------------------------------------------------------------------------
    // Parse and expand input parents

    // if args.parents supplied on the CLI
    if let Some(parents) = &args.parents {
        info!("Parsing input parents: {:?}", &parents);
        let parents_expanded = dataset.expand_populations(parents)?;
        parent_search_populations.retain(|pop| parents_expanded.contains(pop));
        args.parents = Some(parents_expanded);
    }

    // ------------------------------------------------------------------------
    // Knockout

    // if args.knockout supplied on the CLI
    if let Some(knockout) = &args.knockout {
        info!("Performing dataset knockout: {knockout:?}");

        // Check to make sure wildcards are used in all knockouts
        // Weird things happend if you don't!
        let knockout_no_wildcard =
            knockout.iter().filter(|p| !p.contains('*')).collect_vec();
        if !knockout_no_wildcard.is_empty() {
            warn!("Proceed with caution! Weird things can happen when you run a knockout without descendants (no '*'): {knockout_no_wildcard:?}.")
        }

        // Expanded populations (in case wildcard * is provided)
        let knockout_expanded = dataset.expand_populations(knockout)?;
        debug!("Expanded dataset knockout: {knockout:?}");
        debug!("Removing knockout populations from the fasta.");
        dataset.populations.retain(|id, _| !knockout_expanded.contains(id));

        debug!("Removing knockout populations from the mutations.");
        dataset.mutations.iter_mut().for_each(|(_sub, populations)| {
            populations.retain(|p| !knockout_expanded.contains(p));
        });

        if !dataset.phylogeny.is_empty() {
            for p in &knockout_expanded {
                dataset.phylogeny.remove(p)?;
            }
        }

        parent_search_populations.retain(|pop| !knockout_expanded.contains(pop));
        args.knockout = Some(knockout_expanded);
    }

    // ------------------------------------------------------------------------
    // Input Parsing
    // ------------------------------------------------------------------------

    // To process --populations and --alignment in parallel simultaneously
    // we chain together fasta::Records.

    // ------------------------------------------------------------------------
    // Populations (--populations)

    let mut num_population_records = 0;

    // get sequences for each request population
    let sequences = if let Some(populations) = &args.input.populations {
        info!("Parsing input populations: {populations:?}");
        dataset
            .expand_populations(populations)?
            .into_iter()
            .map(|p| {
                let sequence = dataset.populations.get(&p).unwrap();
                let id = format!("population_{}", sequence.id);
                num_population_records += 1;
                format!(">{id}\n{}", sequence.seq.iter().join(""))
            })
            .join("\n")
    } else {
        String::new()
    };

    let mut reader = fasta::Reader::new(sequences.as_bytes());
    let population_records = reader.records();

    // ------------------------------------------------------------------------
    // Alignment (--alignment)

    let temp_file = tempfile::NamedTempFile::new()?;
    let path = match &args.input.alignment {
        Some(alignment) => alignment,
        None => temp_file.path(),
    };

    // read alignment to count records
    let mut num_alignment_records = 0;
    if let Some(alignment) = &args.input.alignment {
        info!("Parsing input alignment: {alignment:?}");

        // count number of records in alignment
        let progress_bar_style = ProgressStyle::with_template(
            "{bar:40} {pos}/{len} ({percent}%) | Sequences / Second: {per_sec} | Elapsed: {elapsed_precise}"
        ).wrap_err("Failed to create progress bar from template.")?;
        let progress_bar = ProgressBar::new(0_u64);
        progress_bar.set_style(progress_bar_style);

        let mut reader = File::open(path).map(BufReader::new).map(fasta::Reader::new)?;
        num_alignment_records = reader
            .records()
            .par_bridge()
            .inspect(|_| {
                progress_bar.inc_length(1);
                progress_bar.inc(1);
            })
            .count();
        progress_bar.finish();
    }

    // read alignment into interator
    let mut reader = File::open(path).map(BufReader::new).map(fasta::Reader::new)?;
    let alignment_records = reader.records();

    // ------------------------------------------------------------------------
    // Combine

    let records = population_records.chain(alignment_records);
    let num_records = num_population_records + num_alignment_records;

    // ------------------------------------------------------------------------
    // Linelist
    // ------------------------------------------------------------------------

    let linelist_path = args.output_dir.join("linelist.tsv");
    info!("Initializing linelist: {linelist_path:?}");

    // get the delimiter based on the extension (ex. tsv => tab)
    let linelist_delim = utils::path_to_delim(&linelist_path)?.to_string();

    // create linelist output file and write headers
    let headers = export::linelist_headers();
    let line = format!("{}\n", headers.join(&linelist_delim));
    let mut file = File::create(&linelist_path)?;
    file.write_all(line.as_bytes())?;

    // append to linelist progressively, use Mutex lock for parallel access
    let linelist_file = OpenOptions::new().append(true).open(linelist_path)?;
    let linelist_file = Mutex::new(linelist_file);

    // ------------------------------------------------------------------------
    // Detect Recombination
    // ------------------------------------------------------------------------

    info!("Detecting recombination.");

    // progress bar style
    let progress_bar_style = ProgressStyle::with_template(
        "{bar:40} {pos}/{len} ({percent}%) | Sequences / Second: {per_sec} | Elapsed: {elapsed_precise} | ETA: {eta_precise}"
    ).wrap_err("Failed to create progress bar from template.")?;

    let progress_bar = ProgressBar::new(num_records as u64);
    progress_bar.set_style(progress_bar_style);

    // if an error is encountered, it will stop
    let results: Vec<_> = records
        .par_bridge()
        .filter_map(|result| result.ok())
        .map(|r| {
            //progress_bar.inc_length(1);
            let seq = Sequence::from_record(r, Some(&dataset.reference), &args.mask)?;
            let (best_match, recombination) =
                search(&seq, &dataset, &parent_search_populations, args)?;

            // append to linelist
            let linelist_table = export::linelist(&best_match, &recombination, &dataset)?;
            let row = &linelist_table.rows[0];
            let line = format!("{}\n", row.join(&linelist_delim));
            linelist_file.lock().unwrap().write_all(line.as_bytes())?;

            progress_bar.inc(1);

            let result = (seq, best_match, recombination);
            Ok(result)
        })
        .collect::<Result<_, Report>>()?;

    progress_bar.finish();

    // ------------------------------------------------------------------------
    // Export Barcodes (multiple, collected by recombinant)

    let outdir_barcodes = args.output_dir.join("barcodes");
    info!("Exporting recombination barcodes: {outdir_barcodes:?}");
    create_dir_all(&outdir_barcodes)?;

    // get unique keys of recombinants identified
    let unique_keys = results
        .iter()
        .filter_map(|(_s, _b, r)| {
            (*r.unique_key != String::new()).then_some(&r.unique_key)
        })
        .unique()
        .collect_vec();

    if unique_keys.is_empty() {
        warn!("No recombination detected, no barcodes will be outputted.");
    }

    // iterate through all unique recombinants
    unique_keys.into_iter().try_for_each(|unique_key| {
        // filter to Recombinations and Sequences that match this unique key
        let (sequences, recombinations): (Vec<_>, Vec<_>) = results
            .iter()
            .filter_map(|(s, _b, r)| (r.unique_key == *unique_key).then_some((s, r)))
            .unzip();
        // combine all the sample barcode tables
        let barcode_table = recombination::combine_tables(
            &sequences,
            &recombinations,
            &dataset.reference,
        )?;
        let barcode_table_path = outdir_barcodes.join(format!("{unique_key}.tsv"));
        barcode_table.write(&barcode_table_path)?;
        Ok::<(), Report>(())
    })?;

    info!("Done.");
    Ok(())
}

pub fn search(
    sequence: &Sequence,
    dataset: &Dataset,
    populations: &[String],
    args: &cli::run::Args,
) -> Result<(SearchResult, Recombination), Report> {
    // initialize with default results, regardless of whether our
    // searches "succeed", we're going to return standardized data
    // structures to build our exports upon (ex. linelist columns)
    // which will include the "negative" results
    let mut best_match = SearchResult::new(sequence);
    let mut recombination = Recombination::new();

    // ------------------------------------------------------------------------
    // Best Match (Consensus)
    //
    // search for the best match in the dataset to this sequence.
    // this will represent the consensus population call.

    debug!(
        "Identifying best match (consensus population): {}",
        &sequence.id
    );
    let search_result = dataset.search(sequence, None, None);

    // if we found a match, proceed with recombinant search
    if let Ok(search_result) = search_result {
        best_match = search_result;

        debug!("Searching for recombination parents: {}", &sequence.id);
        let parent_search = recombination::search::all_parents(
            sequence,
            dataset,
            &mut best_match,
            populations,
            args,
        );
        match parent_search {
            Ok(search_result) => recombination = search_result,
            Err(e) => debug!("Parent search did not succeed. {e}"),
        }
    }
    // what to do if not a single population matched?
    // I think this can only happen if the sample is essentially the reference genome
    else {
        // temporary handling for sars-cov-2 root population B
        if dataset.name == Name::SarsCov2 && sequence.id == "population_B" {
            best_match.consensus_population = "B".to_string();
        } else {
            warn!("No dataset matches were found: {}", &sequence.id);
        }
    }

    Ok((best_match, recombination))
}
