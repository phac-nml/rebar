use rebar::cli;
use rebar::dataset::attributes::{Name, Tag};
use rebar::dataset::download;
use rebar::plot::plot;
use rebar::run::run;

use color_eyre::eyre::{Report, Result};
use std::path::PathBuf;
use std::str::FromStr;

#[tokio::test]
async fn toy1() -> Result<(), Report> {
    let output_dir = PathBuf::from("output").join("tests").join("toy1");

    // Dataset Download
    let mut args = cli::dataset::download::Args {
        name: Name::Toy1,
        tag: Tag::from_str("custom")?,
        output_dir: output_dir.join("dataset"),
        summary: None,
    };
    download::dataset(&mut args).await?;

    // Run
    let mut args = cli::run::Args {
        population: Some("*".to_string()),
        dataset_dir: output_dir.join("dataset"),
        output_dir: output_dir.join("run"),
        mask: vec![0, 0],
        min_length: 3,
        ..Default::default()
    };
    run(&mut args)?;

    // Plot
    let args = cli::plot::Args {
        dataset_dir: output_dir.join("dataset"),
        run_dir: output_dir.join("run"),
        ..Default::default()
    };
    plot(&args)?;

    Ok(())
}

#[tokio::test]
async fn sarscov2_populations() -> Result<(), Report> {
    let output_dir =
        PathBuf::from("output").join("tests").join("sarscov2").join("populations");

    // Dataset Download
    let mut args = cli::dataset::download::Args {
        name: Name::SarsCov2,
        tag: Tag::from_str("2023-11-17")?,
        output_dir: output_dir.join("dataset"),
        summary: None,
    };
    download::dataset(&mut args).await?;

    // Run
    let mut args = cli::run::Args {
        population: Some("AY.4.2*,BA.5.2,XBC.1.6*,XBB.1.5.1,XBL".to_string()),
        dataset_dir: output_dir.join("dataset"),
        output_dir: output_dir.join("run"),
        ..Default::default()
    };
    run(&mut args)?;

    // Plot
    let args = cli::plot::Args {
        dataset_dir: output_dir.join("dataset"),
        run_dir: output_dir.join("run"),
        ..Default::default()
    };
    plot(&args)?;

    Ok(())
}
