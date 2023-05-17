use crate::sequence::Sequence;
use crate::traits::ToYaml;
use color_eyre::eyre::{eyre, Report, WrapErr};
use color_eyre::section::Section;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::create_dir_all;
use std::io::Read;
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

        create_dir_all(output_dir)?;
        dataset.download_reference(output_dir).await?;
        dataset.download_populations(output_dir).await?;
        Ok(dataset)
    }

    pub async fn download_reference(
        &self,
        output_dir: &Path,
    ) -> Result<(), Report> {
        let url = format!(
            "{}/{}/references/{}/versions/{}/files/reference.fasta",
            NEXTCLADE_DATA_URL, self.name, self.reference, self.tag,
        );
        let output_path = output_dir.join("reference.fasta");

        download_file(&url, &output_path).await?;

        Ok(())
    }

    pub async fn download_populations(
        &self,
        output_dir: &Path,
    ) -> Result<(), Report> {
        println!("download_populations");
        let url = match self.name {
            Name::SarsCov2 => "https://raw.githubusercontent.com/corneliusroemer/pango-sequences/main/data/pango-consensus-sequences_genome-nuc.fasta.zst",
            Name::RsvA => return Err(eyre!("Not implemented yet")),
            Name::RsvB => return Err(eyre!("Not implemented yet")),
        };

        let file_ext = Path::new(url).extension().unwrap();

        let output_path = match file_ext.to_str().unwrap() {
            "zst" => output_dir.join("populations.fasta.zst"),
            _ => output_dir.join("populations.fasta"),
        };
        download_file(url, &output_path).await?;

        let input_path = output_path;
        let output_path = output_dir.join("populations.fasta");
        decompress_file(&input_path, &output_path)?;

        Ok(())
    }
}

pub async fn download_file(
    url: &str,
    output_path: &PathBuf,
) -> Result<(), Report> {
    let response = reqwest::get(url).await?;
    if response.status() != 200 {
        return Err(eyre!(
            "Unable to download file: {url}\nStatus code {}.",
            response.status()
        ));
    }
    let content = response.text().await?;

    std::fs::write(output_path, content)
        .wrap_err(format!("Unable to write file: {:?}", output_path))
}

pub fn decompress_file(
    input: &PathBuf,
    output: &PathBuf,
) -> Result<(), Report> {
    let ext = input.extension().unwrap();

    match ext.to_str().unwrap() {
        "zst" => {
            let mut f = std::fs::File::open(input)?;
            let mut buffer = String::new();
            f.read_to_string(&mut buffer)?;
            std::fs::write(output, buffer)
                .wrap_err(format!("Unable to write file: {:?}", output))?;
            std::fs::remove_file(input)?;
        }
        _ => println!("Not implemented yet"),
    };

    Ok(())
}
