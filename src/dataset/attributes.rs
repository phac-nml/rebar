use crate::utils;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use color_eyre::Help;
use itertools::Itertools;
use log::info;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::str::FromStr;

// ----------------------------------------------------------------------------
// Dataset Name

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize, PartialEq)]
pub enum Name {
    #[default]
    SarsCov2,
    RsvA,
    RsvB,
    Unknown,
}

impl fmt::Display for Name {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Name::SarsCov2 => String::from("sars-cov-2"),
            Name::RsvA => String::from("rsv-a"),
            Name::RsvB => String::from("rsv-b"),
            Name::Unknown => String::from("unknown"),
        };

        write!(f, "{}", name)
    }
}

impl FromStr for Name {
    type Err = Report;

    fn from_str(name: &str) -> Result<Self, Report> {
        let name = match name {
            "sars-cov-2" => Name::SarsCov2,
            "rsv-a" => Name::RsvA,
            "rsv-b" => Name::RsvB,
            "unknown" => Name::Unknown,
            _ => Err(eyre!("Unknown dataset name: {name}"))
                .suggestion("Please choose from:")?,
        };

        Ok(name)
    }
}

// ----------------------------------------------------------------------------
// Dataset Tag

#[derive(Clone, Debug, Default, Deserialize, Serialize, PartialEq)]
pub enum Tag {
    #[default]
    Latest,
    Nightly,
    Archive(String),
    Unknown,
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let tag = match self {
            Tag::Latest => String::from("latest"),
            Tag::Nightly => String::from("nightly"),
            Tag::Archive(tag) => tag.to_owned(),
            Tag::Unknown => String::from("unknown"),
        };

        write!(f, "{}", tag)
    }
}

impl FromStr for Tag {
    type Err = Report;

    fn from_str(tag: &str) -> Result<Tag, Report> {
        let tag = match tag {
            "latest" => Tag::Latest,
            "nightly" => Tag::Nightly,
            "unknown" => Tag::Unknown,
            _ => Tag::Archive(String::from(tag)),
        };

        Ok(tag)
    }
}

// ----------------------------------------------------------------------------
// Dataset Summary

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Summary {
    pub tag: Tag,
    pub name: Name,
}

impl Summary {
    /// import summary from specified format
    pub fn import(
        dataset_dir: &Path,
        format: SummaryImportFormat,
    ) -> Result<Summary, Report> {
        //  import path
        let mut import_path = dataset_dir.join("summary");
        import_path.set_extension(format.extension());

        let summary: Summary = match format {
            SummaryImportFormat::Json => {
                let summary = std::fs::read_to_string(import_path)
                    .wrap_err_with(|| "Couldn't read summary {import_path:?}.")?;
                serde_json::from_str(&summary)?
            }
        };

        Ok(summary)
    }

    /// Export summary to specified format
    pub fn export(
        &self,
        dataset_dir: &Path,
        format: SummaryExportFormat,
    ) -> Result<(), Report> {
        // output path
        let mut output_path = dataset_dir.join("summary");
        output_path.set_extension(format.extension());
        info!("Exporting summary to {format}: {output_path:?}");

        // format conversion
        let output = match format {
            SummaryExportFormat::Json => serde_json::to_string_pretty(self)
                .wrap_err_with(|| format!("Failed to export summary to {format}"))?,
        };

        // Write to file
        let mut file = File::create(&output_path).wrap_err_with(|| {
            format!("Failed to access output summary path {output_path:?}")
        })?;
        file.write_all(output.as_bytes())
            .wrap_err_with(|| format!("Failed to write summary to {output_path:?}"))?;

        Ok(())
    }
}

// ----------------------------------------------------------------------------
// Dataset Summary Export Format

pub enum SummaryExportFormat {
    Json,
}

impl SummaryExportFormat {
    /// get extension from format
    pub fn extension(&self) -> String {
        match self {
            SummaryExportFormat::Json => String::from("json"),
        }
    }
}

impl std::fmt::Display for SummaryExportFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SummaryExportFormat::Json => write!(f, "json"),
        }
    }
}

// ----------------------------------------------------------------------------
// Dataset Summary Import Format

pub enum SummaryImportFormat {
    Json,
}

impl SummaryImportFormat {
    /// get extension from format
    pub fn extension(&self) -> String {
        match self {
            SummaryImportFormat::Json => String::from("json"),
        }
    }
}

impl std::fmt::Display for SummaryImportFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SummaryImportFormat::Json => write!(f, "json"),
        }
    }
}

// ----------------------------------------------------------------------------
// Dataset Annotations

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Annotations {
    pub gene: Vec<String>,
    pub abbreviation: Vec<String>,
    pub start: Vec<usize>,
    pub end: Vec<usize>,
}

impl Annotations {
    pub fn to_table(&self) -> Result<utils::table::Table, Report> {
        let mut table = utils::table::Table::new();

        // headers
        table.headers = vec!["gene", "abbreviation", "start", "end"]
            .into_iter()
            .map(String::from)
            .collect_vec();

        // rows
        let mut rows = Vec::new();
        for i in 0..(self.gene.len()) {
            let row = vec![
                self.gene[i].clone(),
                self.abbreviation[i].clone(),
                self.start[i].to_string(),
                self.end[i].to_string(),
            ];
            rows.push(row);
        }

        table.rows = rows;

        Ok(table)
    }
}
