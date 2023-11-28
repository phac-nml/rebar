use crate::cli;
use crate::dataset::attributes::{Name, Tag};
use crate::dataset::download;
use color_eyre::eyre::{Report, Result};
use std::path::PathBuf;
use std::str::FromStr;

#[tokio::test]
async fn example_1() -> Result<(), Report> {
    // ------------------------------------------------------------------------
    // Dataset Download

    let name = Name::SarsCov2;
    let tag = Tag::from_str("2023-11-17")?;
    let output_dir = PathBuf::from("test").join("example_1").join("dataset");
    let summary = None;

    let mut args = cli::dataset::download::Args {
        name,
        tag,
        output_dir: output_dir.clone(),
        summary,
    };

    download::dataset(&mut args).await?;

    // // ------------------------------------------------------------------------
    // // Run

    // let population = Some("AY.4.2*,BA.5.2,XBC.1.6*,XBB.1.5.1,XBL".to_string());
    // let dataset_dir = output_dir;
    // let output_dir = PathBuf::from("test").join("example_1").join("run");

    // let mut args = cli::run::Args {
    //     population,
    //     dataset_dir: dataset_dir.clone(),
    //     output_dir: output_dir.clone(),
    //     ..Default::default()
    // };

    // run(&mut args)?;

    // ------------------------------------------------------------------------
    // Plot

    Ok(())
}
