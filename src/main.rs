use clap::Parser;
use color_eyre::eyre::{Report, Result};
use log::{debug, info};
use rebar::cli::verbosity::Verbosity;
use rebar::cli::{Cli, Command};
use rebar::dataset::Dataset;
use rebar::query::Query;
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
            Dataset::download(&name, &tag, &output_dir).await?;
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
            // Load dataset
            let dataset = Dataset::load(&dataset_dir, mask)?;
            // Load the query alignment
            let query = Query::load(alignment, &dataset, mask)?;
            info!("Identifying consensus and parent populations.");
            for (id, sequence) in query.sequences {
                //if id != "XBB.1.16" && id != "BM.1.1.1" {
                if id != "XBB.1.16" {
                    continue;
                }

                // best match ie. consensus population
                let exclude_populations = None;
                let include_populations = None;
                debug!("sequence: {id}");
                let best_match =
                    dataset.find_best_match(&sequence, exclude_populations, include_populations)?;
                let _parents =
                    dataset.find_parents(sequence, &best_match, max_parents, max_iter, min_consecutive, min_length, min_subs)?;

                // for (i, match_summary) in parents.iter().enumerate() {
                //     debug!("parent_{}: {}", i+1, match_summary.consensus_population);
                //     debug!("{}", match_summary.to_yaml().replace('\n', format!("\n{}", " ".repeat(31)).as_str()));
                // }
            }
        }
    }

    Ok(())
}
