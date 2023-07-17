use clap::Parser; // CLI argument parsing
use color_eyre::eyre::{Report, Result}; // colorized logging and error handling
use rebar::cli::{Cli, Command, DatasetCommand};

#[tokio::main]
async fn main() -> Result<(), Report> {
    // Parse CLI parameters
    let args = Cli::parse();

    // ------------------------------------------------------------------------
    // CLI Setup

    // initialize color_eyre crate for colorized logs
    color_eyre::install()?;

    // Set logging/verbosity level via RUST_LOG
    std::env::set_var("RUST_LOG", args.verbosity.to_string());

    // initialize env_logger crate for logging/verbosity level
    env_logger::init();

    // check which CLI command we're running (dataset, run, plot)
    match args.command {
        // --------------------------------------------------------------------
        // Dataset
        Command::Dataset(dataset_args) => match dataset_args.command {
            DatasetCommand::List(_) => todo!(),
            DatasetCommand::Download(args) => rebar::dataset_download(args).await?,
        },
        // --------------------------------------------------------------------
        // Run
        Command::Run(args) => rebar::run(*args)?,
        // --------------------------------------------------------------------
        // Plot
        Command::Plot(args) => rebar::plot(*args)?,
    }

    Ok(())
}
