use crate::utils::remote_file::RemoteFile;
use chrono::prelude::*;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use color_eyre::Help;
use indoc::formatdoc;
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
    #[strum(props(implemented = "true"))]
    SarsCov2,
    #[serde(rename = "rsv-a")]
    #[strum(props(implemented = "false"))]
    RsvA,
    #[serde(rename = "rsv-b")]
    #[strum(props(implemented = "false"))]
    RsvB,
    #[default]
    #[serde(rename = "custom")]
    #[strum(props(implemented = "false"))]
    Custom,
}

impl Name {
    pub fn compatibility(&self) -> Result<Compatibility, Report> {
        let mut compatibility = Compatibility::new();
        #[allow(clippy::single_match)]
        match self {
            Name::SarsCov2 => {
                compatibility.dataset.min_date =
                    Some(NaiveDate::parse_from_str("2023-02-09", "%Y-%m-%d")?);
            }
            _ => compatibility.cli.version = Some(">=1.0.0".to_string()),
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
            Name::Custom => String::from("custom"),
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
            "custom" => Name::Custom,
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
                let tag_date = NaiveDate::parse_from_str(tag, "%Y-%m-%d")
                    .wrap_err_with(|| eyre!("Archive tag date is invalid: {tag:?}. Example of a valid Archive tag: 2023-08-17"))?;
                // is it in the future?
                let today = Local::now().date_naive();
                if tag_date > today {
                    return Err(eyre!("Archive tag date is in the future: {tag:?}. Please pick a date on or before today: {today:?}"))?;
                }
                Tag::Archive(tag.to_string())
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
        let tag_date = NaiveDate::parse_from_str(&tag.to_string(), "%Y-%m-%d")?;

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
    pub min_date: Option<NaiveDate>,
    pub max_date: Option<NaiveDate>,
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
            tag: Tag::Custom,
            name: Name::Custom,
            reference: RemoteFile::new(),
            populations: RemoteFile::new(),
            misc: BTreeMap::new(),
        }
    }
    /// Read summary from file.
    pub fn read(path: &Path) -> Result<Summary, Report> {
        let summary = std::fs::read_to_string(path)
            .wrap_err_with(|| format!("Failed to read file: {path:?}."))?;
        let summary = serde_json::from_str(&summary)
            .wrap_err_with(|| format!("Failed to parse file: {path:?}"))?;

        Ok(summary)
    }

    /// Write summary to file.
    pub fn write(&self, path: &Path) -> Result<(), Report> {
        // create output file
        let mut file = File::create(path)
            .wrap_err_with(|| format!("Failed to create file: {path:?}"))?;

        // parse to string
        let output = serde_json::to_string_pretty(self)
            .wrap_err_with(|| format!("Failed to parse: {self:?}"))?;

        // write to file
        file.write_all(format!("{}\n", output).as_bytes())
            .wrap_err_with(|| format!("Failed to write file: {path:?}"))?;

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
