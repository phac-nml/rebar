use crate::utils::remote_file::RemoteFile;
use chrono::prelude::*;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use color_eyre::Help;
use indoc::formatdoc;
use log::debug;
use semver::{Version, VersionReq};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::default::Default;
use std::fmt;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::str::FromStr;
use strum::{EnumIter, EnumProperty};

// ----------------------------------------------------------------------------
// Dataset Name

#[derive(
    Clone, Copy, Debug, Default, EnumIter, EnumProperty, Serialize, Deserialize, PartialEq,
)]
pub enum Name {
    #[serde(rename = "sars-cov-2")]
    #[strum(props(list = "true"))]
    SarsCov2,
    #[serde(rename = "rsv-a")]
    #[strum(props(list = "false"))]
    RsvA,
    #[serde(rename = "rsv-b")]
    #[strum(props(ilist = "false"))]
    RsvB,
    #[default]
    #[strum(props(list = "false"))]
    Unknown,
}

impl Name {
    pub fn compatibility(&self) -> Result<Compatibility, Report> {
        let mut compatibility = Compatibility::new();
        #[allow(clippy::single_match)]
        match self {
            Name::SarsCov2 => {
                compatibility.dataset.min_date =
                    Some(DateTime::parse_from_rfc3339("2023-02-09T12:00:00Z")?.into());
            }
            _ => (),
        }
        Ok(compatibility)
    }
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
    Latest,
    Archive(String),
    #[default]
    Unknown,
    Custom,
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let tag = match self {
            Tag::Latest => String::from("latest"),
            Tag::Archive(tag) => tag.to_owned(),
            Tag::Custom => String::from("custom"),
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
            "unknown" => Tag::Unknown,
            "custom" => Tag::Custom,
            _ => {
                // check if it's an archival date string
                let tag_date = DateTime::parse_from_rfc3339(tag)
                    .wrap_err_with(|| eyre!("Archive tag date has invalid format: {tag:?}. Example of a valid Archive tag: 2023-08-17T12:00:00Z"))?;
                let tag_utc: DateTime<Utc> = tag_date.into();
                let tag_reformat = tag_utc.to_rfc3339_opts(SecondsFormat::Secs, true);
                Tag::Archive(tag_reformat)
            }
        };

        Ok(tag)
    }
}

// ----------------------------------------------------------------------------
// Dataset Compatibility

pub fn check_compatibility(name: &Name, tag: &Tag) -> Result<(), Report> {
    let compatibility = name.compatibility()?;

    // Check CLI Version
    if let Some(cli_version) = compatibility.cli.version {
        let current_version = Version::parse(env!("CARGO_PKG_VERSION"))?;
        let required_version = VersionReq::parse(&cli_version)?;
        if !required_version.matches(&current_version) {
            return Err(eyre!(formatdoc!(
                "CLI version incompatibility.
                Current version {current_version} does not satisfy the {name} dataset requirement {required_version}",
                current_version=current_version.to_string()
                )));
        }
    }
    // Check Tag Dates
    if matches!(tag, Tag::Archive(_)) {
        let tag_date: DateTime<Utc> =
            DateTime::parse_from_rfc3339(&tag.to_string())?.into();

        // Minimum Date
        if let Some(min_date) = compatibility.dataset.min_date {
            if tag_date < min_date {
                return Err(eyre!(formatdoc!(
                    "Date incompatibility.
                    Tag {tag_date:?} does not satisfy the {name} dataset minimum date {min_date:?}"
                )));
            }
        }
        // Maximum Date
        if let Some(max_date) = compatibility.dataset.max_date {
            if tag_date > max_date {
                return Err(eyre!(formatdoc!(
                    "Date incompatibility.
                    Tag {tag_date:?} does not satisfy the {name} dataset maximum date {max_date:?}"
                )));
            }
        }
    }

    Ok(())
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Compatibility {
    pub dataset: DateCompatibility,
    pub cli: CliCompatibility,
}

impl Default for Compatibility {
    fn default() -> Self {
        Self::new()
    }
}

impl Compatibility {
    pub fn new() -> Self {
        Compatibility {
            dataset: DateCompatibility::new(),
            cli: CliCompatibility::new(),
        }
    }
}

// ----------------------------------------------------------------------------
// Date Compatibility

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DateCompatibility {
    pub min_date: Option<DateTime<Utc>>,
    pub max_date: Option<DateTime<Utc>>,
}

impl Default for DateCompatibility {
    fn default() -> Self {
        Self::new()
    }
}

impl DateCompatibility {
    pub fn new() -> Self {
        DateCompatibility {
            min_date: None,
            max_date: None,
        }
    }
}

// ----------------------------------------------------------------------------
// CLI Compatibility

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CliCompatibility {
    pub version: Option<String>,
}

impl Default for CliCompatibility {
    fn default() -> Self {
        Self::new()
    }
}

impl CliCompatibility {
    pub fn new() -> Self {
        CliCompatibility {
            version: Some(">=0.1.0".to_string()),
        }
    }
}

// ----------------------------------------------------------------------------
// Dataset Summary

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Summary {
    pub version: String,
    pub tag: Tag,
    pub name: Name,
    pub reference: RemoteFile,
    pub populations: RemoteFile,
    pub misc: BTreeMap<String, RemoteFile>,
}

impl Default for Summary {
    fn default() -> Self {
        Self::new()
    }
}
impl Summary {
    pub fn new() -> Self {
        Summary {
            version: format!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION")),
            tag: Tag::Unknown,
            name: Name::Unknown,
            reference: RemoteFile::new(),
            populations: RemoteFile::new(),
            misc: BTreeMap::new(),
        }
    }
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
        debug!("Exporting summary to {format}: {output_path:?}");

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
