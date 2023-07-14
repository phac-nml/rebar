use clap::Parser; // CLI argument parsing
use color_eyre::eyre::{Report, Result}; // colorized logging and error handling

#[tokio::main]
async fn main() -> Result<(), Report> {
    // Parse CLI parameters
    let args = rebar::cli::Cli::parse();

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
        rebar::cli::Command::Dataset(dataset_args) => match dataset_args.command {
            rebar::cli::DatasetCommand::List(_dataset_list_args) => todo!(),
            rebar::cli::DatasetCommand::Download(dataset_download_args) => {
                rebar::dataset_download(dataset_download_args).await?
            }
        },
        // --------------------------------------------------------------------
        // Run
        rebar::cli::Command::Run(run_args) => rebar::run(*run_args)?,
        // --------------------------------------------------------------------
        // Plot
        rebar::cli::Command::Plot(plot_args) => rebar::plot(*plot_args)?,
    }

    Ok(())
}
