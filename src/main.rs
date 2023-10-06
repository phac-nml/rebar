use clap::Parser;
use color_eyre::eyre::{Report, Result};
use rebar::cli::{dataset, Cli, Command};

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
        Command::Dataset(args) => match args.command {
            dataset::Command::List(args) => rebar::list_datasets(&args).await?,
            dataset::Command::Download(mut args) => {
                rebar::download_dataset(&mut args).await?
            }
        },
        // --------------------------------------------------------------------
        // Run
        Command::Run(mut args) => rebar::run(&mut args)?,
        // --------------------------------------------------------------------
        // Plot
        Command::Plot(args) => rebar::plot(*args)?,
        // --------------------------------------------------------------------
        // Simulate
        Command::Simulate(args) => rebar::simulate(&args)?,
    }

    Ok(())
}
