use bio::io::fasta;
use clap::Parser;
use color_eyre::eyre::{Report, Result, WrapErr};
use log::{debug, info};
use rebar::cli::{Cli, Command, Verbosity};
use rebar::dataset;
use rebar::sequence::Sequence;
use std::env;
use std::str::FromStr;

fn setup(args: &Cli) -> Result<(), Report> {
    // initialize color_eyre for colorized logs
    color_eyre::install()?;

    let verbosity = match &args.command {
        Command::Dataset { verbosity, .. } => Verbosity::from_str(verbosity)?,
        Command::Run { verbosity, .. } => Verbosity::from_str(verbosity)?,
    };
    // Set default logging level via RUST_LOG
    env::set_var("RUST_LOG", verbosity.to_string());

    // initialize env_logger for different verbosity levels
    env_logger::init();

    Ok(())
}

#[tokio::main]
async fn main() -> Result<()> {
    // Parse CLI parameters
    let args = Cli::parse();

    // Misc setup actions like logging
    setup(&args)?;

    match args.command {
        // Download a dataset
        Command::Dataset {
            name,
            tag,
            output_dir,
            ..
        } => {
            dataset::download(&name, &tag, &output_dir).await?;
        }
        // Run on input alignment
        Command::Run {
            alignment,
            dataset_dir,
            mask,
            max_iter,
            max_parents,
            min_consecutive,
            min_length,
            min_subs,
            ..
        } => {
            // Collect files in dataset_dir into a dataset object
            // This mainly includes parent populations sequences
            //   and optionally a phylogenetic representation.
            let dataset = dataset::load(&dataset_dir, mask)?;

            info!("Loading alignment: {:?}", alignment);
            let alignment_reader =
                fasta::Reader::from_file(&alignment).expect("Unable to read alignment");

            for result in alignment_reader.records() {
                let record = result.wrap_err(format!(
                    "Unable to parse alignment: {:?}",
                    alignment.to_str().unwrap()
                ))?;
                let sequence =
                    Sequence::from_record(record, Some(&dataset.reference), mask)?;

                if sequence.id != "XBB.1.16" {
                    continue;
                }

                // best match ie. consensus population
                let exclude_populations = None;
                let include_populations = None;
                debug!("sequence: {}", sequence.id);
                let best_match = dataset::find_best_match(
                    &dataset,
                    &sequence,
                    exclude_populations,
                    include_populations,
                )?;
                let _parents = dataset::find_parents(
                    &dataset,
                    sequence,
                    &best_match,
                    max_parents,
                    max_iter,
                    min_consecutive,
                    min_length,
                    min_subs,
                )?;
            }
        }
    }

    Ok(())
}
