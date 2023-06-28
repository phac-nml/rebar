use clap::Parser; // CLI argument parsing
use color_eyre::eyre::Result; // colorized logging and error handling

#[tokio::main]
async fn main() -> Result<()> {
    // Parse CLI parameters
    let args = rebar::cli::Cli::parse();

    // --------------------------------------------------------------------
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
        rebar::cli::Command::Dataset(args) => {
            match args.command {
                // list
                rebar::cli::DatasetCommand::List(_) => todo!(),
                // download
                rebar::cli::DatasetCommand::Download(args) => {
                    rebar::dataset::download(&args.name, &args.tag, &args.output_dir)
                        .await?;
                }
            }
        }

        // --------------------------------------------------------------------
        // Run
        rebar::cli::Command::Run(args) => {
            rebar::run(*args)?;
        }
        // --------------------------------------------------------------------
        // Plot
        rebar::cli::Command::Plot(args) => {
            rebar::plot(*args)?;
        }
    }

    Ok(())
}
