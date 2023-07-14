use crate::cli;
use crate::dataset;
use color_eyre::eyre::{Report, Result, WrapErr};
use log::{debug, info};
use serde::{Deserialize, Serialize};
use std::default::Default;
use std::fs::File;
use std::io::Write;
use std::path::Path;

// ----------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// EdgeCase

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct EdgeCase {
    population: String,
    args: cli::RunArgs,
}

impl Default for EdgeCase {
    fn default() -> Self {
        Self::new()
    }
}

impl EdgeCase {
    pub fn new() -> Self {
        EdgeCase {
            population: String::new(),
            args: cli::RunArgs::new(),
        }
    }
}

// ----------------------------------------------------------------------------
// EdgeCaseExportFormat

pub enum EdgeCaseExportFormat {
    Json,
}

impl EdgeCaseExportFormat {
    pub fn extension(&self) -> String {
        match self {
            EdgeCaseExportFormat::Json => String::from("json"),
        }
    }
}

impl std::fmt::Display for EdgeCaseExportFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            EdgeCaseExportFormat::Json => write!(f, "json"),
        }
    }
}

// ----------------------------------------------------------------------------
// EdgeCaseImportFormat

pub enum EdgeCaseImportFormat {
    Json,
}

impl EdgeCaseImportFormat {
    pub fn extension(&self) -> String {
        match self {
            EdgeCaseImportFormat::Json => String::from("json"),
        }
    }
}

impl std::fmt::Display for EdgeCaseImportFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            EdgeCaseImportFormat::Json => write!(f, "json"),
        }
    }
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

pub fn from_dataset_name(
    name: &dataset::attributes::Name,
) -> Result<Vec<EdgeCase>, Report> {
    let edge_case: Vec<EdgeCase> = match name {
        dataset::attributes::Name::SarsCov2 => create_sarscov2_edge_cases()?,
        _ => Vec::new(),
    };

    Ok(edge_case)
}

pub fn create_sarscov2_edge_cases() -> Result<Vec<EdgeCase>, Report> {
    info!("Creating SARS-CoV-2 edge cases.");
    let mut edge_cases: Vec<EdgeCase> = Vec::new();

    // default args to share between all recombinants
    let default_args = cli::RunArgs::default();

    // --------------------------------------------------------------------
    // XCF
    // XCF is XBB and FE.1 (XBB.1.18.1) with no unique subs from XBB

    debug!("Creating XCF edge case.");
    let mut xcf = EdgeCase {
        population: "XCF".to_string(),
        args: default_args,
    };
    xcf.args.min_subs = 0;
    edge_cases.push(xcf);

    Ok(edge_cases)
}

/// import edge cases from specified format
pub fn import(
    dataset_dir: &Path,
    format: EdgeCaseImportFormat,
) -> Result<Vec<EdgeCase>, Report> {
    //  import path
    let mut import_path = dataset_dir.join("edge_cases");
    import_path.set_extension(format.extension());

    let edge_cases: Vec<EdgeCase> = match format {
        EdgeCaseImportFormat::Json => {
            let edge_cases = std::fs::read_to_string(import_path)
                .wrap_err_with(|| "Couldn't read edge cases from {import_path:?}.")?;
            serde_json::from_str(&edge_cases)?
        }
    };

    Ok(edge_cases)
}

/// Export edge cases to specified format
pub fn export(
    edge_cases: &[EdgeCase],
    dataset_dir: &Path,
    format: EdgeCaseExportFormat,
) -> Result<(), Report> {
    // output path
    let mut output_path = dataset_dir.join("edge_cases");
    output_path.set_extension(format.extension());
    info!("Exporting edge cases to {format}: {output_path:?}");

    // format conversion
    let output = match format {
        EdgeCaseExportFormat::Json => serde_json::to_string_pretty(edge_cases)
            .wrap_err_with(|| format!("Failed to export edge cases to {format}"))?,
    };

    // Write to file
    let mut file = File::create(&output_path).wrap_err_with(|| {
        format!("Failed to access output edge cases path {output_path:?}")
    })?;
    file.write_all(output.as_bytes())
        .wrap_err_with(|| format!("Failed to write edge cases to {output_path:?}"))?;

    Ok(())
}
