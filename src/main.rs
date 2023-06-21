use bio::io::fasta;
use clap::Parser;
use color_eyre::eyre::{Report, Result, WrapErr};
use log::{debug, info};
use rebar::cli;
use rebar::dataset;
use rebar::sequence::Sequence;
use std::env;
use std::str::FromStr;

#[tokio::main]
async fn main() -> Result<()> {

    // Parse CLI parameters
    let args = cli::RebarCli::parse();
    // initialize color_eyre for colorized logs
    color_eyre::install()?;
    // Set default logging level via RUST_LOG
    env::set_var("RUST_LOG", args.verbosity.to_string());
    // initialize env_logger for different verbosity levels
    env_logger::init();

    match args.command {
        // --------------------------------------------------------------------
        // dataset

        cli::RebarCommand::Dataset(args) => {

            match args.command {
            
                // list
                cli::RebarDatasetCommand::List(_) => todo!(),
                // download
                cli::RebarDatasetCommand::Download(args) => {
                    dataset::download(&args.name, &args.tag, &args.output_dir).await?;
                },
            }
        },

        // --------------------------------------------------------------------
        // run

        cli::RebarCommand::Run(args) => {
            rebar::run(*args)?;          
        },        
    }

    // match args.command {
    //     // Download a dataset
    //     Command::Dataset {
    //         name,
    //         tag,
    //         output_dir,
    //         ..
    //     } => {
    //         dataset::download(&name, &tag, &output_dir).await?;
    //     }
    //     // Run on input alignment
    //     Command::Run {
    //         input,
    //         dataset_dir,
    //         mask,
    //         max_iter,
    //         max_parents,
    //         min_consecutive,
    //         min_length,
    //         min_subs,
    //         ..
    //     } => {
            // Check if the input is a path to alignment or population name
            // println!("{input:?}");
            // let (alignment, population) = (input.alignment, input.population);

            // // Collect files in dataset_dir into a dataset object
            // // This mainly includes parent populations sequences
            // //   and optionally a phylogenetic representation.
            // let dataset = dataset::load(&dataset_dir, mask)?;

            // info!("Loading alignment: {:?}", alignment);
            // let alignment_reader =
            //     fasta::Reader::from_file(&alignment).expect("Unable to read alignment");

            // for result in alignment_reader.records() {
            //     let record = result.wrap_err(format!(
            //         "Unable to parse alignment: {:?}",
            //         alignment.to_str().unwrap()
            //     ))?;
            //     let sequence =
            //         Sequence::from_record(record, Some(&dataset.reference), mask)?;

            //     if sequence.id != "XBB.1.16" {
            //         continue;
            //     }

            //     // best match ie. consensus population
            //     debug!("sequence: {}", sequence.id);
            //     let best_match = dataset.search( &sequence, None, None)?;
            //     let _parents = dataset::find_parents(
            //         &dataset,
            //         sequence,
            //         &best_match,
            //         max_parents,
            //         max_iter,
            //         min_consecutive,
            //         min_length,
            //         min_subs,
            //     )?;
            // }
    //     }
    // }

    Ok(())
}
