//use std::path::Path;

use crate::sequence::Sequence;
use crate::traits::ToYaml;
use color_eyre::eyre::{eyre, Report};
use color_eyre::section::Section;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::default::Default;

#[derive(Debug, Default, PartialEq, Deserialize, Serialize)]
enum DatasetName {
    #[default]
    RsvA,
    RsvB,
    SarsCov2,
}

impl std::fmt::Display for DatasetName {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let lowercase = format!("{:?}", self).to_lowercase();
        write!(f, "{lowercase}")
    }
}

#[derive(Debug, Default, PartialEq, Deserialize, Serialize)]
enum DatasetVersion {
    #[default]
    Latest,
    Nightly,
    // Archive,
}

impl std::fmt::Display for DatasetVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Convert to lowercase for RUST_LOG env var compatibility
        let lowercase = format!("{:?}", self).to_lowercase();
        write!(f, "{lowercase}")
    }
}

#[derive(Debug, Default, Deserialize, Serialize)]
pub struct Dataset {
    name: DatasetName,
    version: DatasetVersion,
    populations: BTreeMap<String, Sequence>,
}

impl std::fmt::Display for Dataset {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "name: {}, version: {}", self.name, self.version)
    }
}

impl ToYaml for Dataset {}

impl Dataset {
    pub fn new(name: String, version: String) -> Result<Dataset, Report> {
        let dataset_name = match name.as_str() {
            "rsv-a" => DatasetName::RsvA,
            "rsv-b" => DatasetName::RsvB,
            "sars-cov-2" => DatasetName::SarsCov2,
            _ => Err(eyre!("Unknown dataset name: {name}"))
                .suggestion("Please choose from: rsva, rsvb, sars-cov2")?,
        };

        let dataset_version = match version.as_str() {
            "latest" => DatasetVersion::Latest,
            "nightly" => DatasetVersion::Nightly,
            _ => Err(eyre!("Unknown dataset version: {version}"))
                .suggestion("Please choose from: latest, nightly")?,
        };

        let dataset = Dataset {
            name: dataset_name,
            version: dataset_version,
            populations: BTreeMap::new(),
        };

        Ok(dataset)
    }

    // pub fn import_sequences()
}
