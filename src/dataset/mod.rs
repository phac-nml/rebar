use crate::sequence::Sequence;
use crate::traits::ToYaml;
use color_eyre::eyre::{eyre, Report, WrapErr};
use color_eyre::section::Section;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::create_dir_all;
use std::path::{Path, PathBuf};

const NEXTCLADE_DATA_URL: &str = "https://raw.githubusercontent.com/nextstrain/nextclade_data/master/data/datasets";

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
pub enum Name {
    RsvA,
    RsvB,
    SarsCov2,
}

impl std::fmt::Display for Name {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Name::RsvA => write!(f, "rsv-a"),
            Name::RsvB => write!(f, "rsv-b"),
            Name::SarsCov2 => write!(f, "sars-cov-2"),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
pub enum Tag {
    Latest,
    Nightly,
    Archive(String),
}

impl std::fmt::Display for Tag {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Tag::Latest => write!(f, "latest"),
            Tag::Nightly => write!(f, "nightly"),
            Tag::Archive(tag) => write!(f, "{tag}"),
        }
    }
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Dataset {
    name: Name,
    tag: Tag,
    reference: String,
    populations: BTreeMap<String, Sequence>,
}

impl std::fmt::Display for Dataset {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "name: {}, tag: {}", self.name, self.tag)
    }
}

impl ToYaml for Dataset {}

impl Dataset {
    pub async fn new(
        name: &String,
        tag: &String,
        reference: &String,
        output_dir: &PathBuf,
    ) -> Result<Dataset, Report> {
        let name = match name.as_str() {
            "rsv-a" => Name::RsvA,
            "rsv-b" => Name::RsvB,
            "sars-cov-2" => Name::SarsCov2,
            _ => Err(eyre!("Unknown dataset name: {name}"))
                .suggestion("Please choose from: rsv-a, rsv-b, sars-cov2")?,
        };

        let tag = match tag.as_str() {
            "latest" => Tag::Latest,
            "nightly" => Tag::Nightly,
            _ => Tag::Archive(tag.to_string()),
        };

        let dataset = Dataset {
            name,
            tag,
            reference: reference.to_string(),
            populations: BTreeMap::new(),
        };

        dataset.download_reference(output_dir).await?;
        Ok(dataset)
    }

    pub async fn download_reference(
        &self,
        output_dir: &PathBuf,
    ) -> Result<(), Report> {
        let url = format!(
            "{}/{}/references/{}/versions/{}/files/reference.fasta",
            NEXTCLADE_DATA_URL, self.name, self.reference, self.tag,
        );

        download_file(&url, output_dir).await?;

        Ok(())
    }

    // pub async fn download_populations(&self, output_dir: &PathBuf) -> Result<(), Report> {

    //     println!("download_populations");
    //     Ok(())
    // }
}

pub async fn download_file(
    url: &str,
    output_dir: &PathBuf,
) -> Result<(), Report> {
    let file_name = Path::new(url).file_name().unwrap();
    let file_path = output_dir.join(file_name);

    let response = reqwest::get(url).await?;
    if response.status() != 200 {
        return Err(eyre!(
            "Unable to download file: {url}\nStatus code {}.",
            response.status()
        ));
    }
    let content = response.text().await?;

    create_dir_all(output_dir)?;
    std::fs::write(&file_path, content)
        .wrap_err(format!("Unable to write file: {:?}", file_path))
}
