use crate::sequence::Sequence;
use crate::traits::ToYaml;
use color_eyre::eyre::{eyre, Report, WrapErr};
use color_eyre::section::Section;
use log::info;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::create_dir_all;
use std::io::Read;
use std::path::{Path, PathBuf};

// const NEXTCLADE_DATA_URL: &str = "https://raw.githubusercontent.com/nextstrain/nextclade_data/master/data/datasets";
const SARSCOV2_REFERENCE_URL: &str = "https://raw.githubusercontent.com/nextstrain/ncov/master/data/references_sequences.fasta";
const SARSCOV2_POPULATIONS_URL: &str = "https://raw.githubusercontent.com/corneliusroemer/pango-sequences/main/data/pango-consensus-sequences_genome-nuc.fasta.zst";

// At minimum, we need a reference and aligned sequences fasta
// Eventually, we might need to wrangle the phylogeny
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

// Just implement nightly first, we'll figure out latest and archive later
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
pub enum Tag {
    Nightly,
    Latest,
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
    reference: Sequence,
    populations: BTreeMap<String, Sequence>,
}

impl std::fmt::Display for Dataset {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "name: {}, tag: {}", self.name, self.tag)
    }
}

impl ToYaml for Dataset {}

impl Dataset {
    pub async fn download(
        name: &String,
        tag: &String,
        output_dir: &Path,
    ) -> Result<(), Report> {
        let name = match name.as_str() {
            "rsv-a" | "rsv-b" => Err(eyre!("Dataset is not implemented yet: {name}"))?,
            "sars-cov-2" => Name::SarsCov2,
            _ => Err(eyre!("Unknown dataset name: {name}"))
                .suggestion("Please choose from:")?,
        };

        let _tag = match tag.as_str() {
            "latest" => Tag::Latest,
            "nightly" => Tag::Nightly,
            _ => Tag::Archive(tag.to_string()),
        };

        if !output_dir.exists() {
            create_dir_all(output_dir)?;
            info!("Creating output directory: {:?}", output_dir);
        }

        // --------------------------------------------------------------------
        // Download Reference
        // --------------------------------------------------------------------
        let url = match name {
            Name::SarsCov2 => SARSCOV2_REFERENCE_URL.to_string(),
            _ => {
                return Err(eyre!(
                    "Downloading the {name} dataset is not implemented yet."
                ))
            }
        };
        let ext = Path::new(&url).extension().unwrap().to_str().unwrap();
        let output_path = match ext {
            "fasta" | "fa" => output_dir.join("reference.fasta"),
            _ => output_dir.join(format!("reference.fasta.{ext}")),
        };
        info!("Downloading reference: {} to {:?}", url, output_path);
        download_file(&url, &output_path).await?;
        if ext != "fasta" && ext != "fa" {
            let output_path_decompress = output_dir.join("reference.fasta");
            decompress_file(&output_path, &output_path_decompress)?;
        }

        // --------------------------------------------------------------------
        // Download Populations
        // --------------------------------------------------------------------
        let url = match name {
            Name::SarsCov2 => SARSCOV2_POPULATIONS_URL.to_string(),
            _ => {
                return Err(eyre!(
                    "Downloading the {name} dataset is not implemented yet."
                ))
            }
        };
        let ext = Path::new(&url).extension().unwrap().to_str().unwrap();
        let output_path = match ext {
            "fasta" | "fa" => output_dir.join("populations.fasta"),
            _ => output_dir.join(format!("populations.fasta.{ext}")),
        };
        info!("Downloading populations: {} to {:?}", url, output_path);
        download_file(&url, &output_path).await?;
        if ext != "fasta" && ext != "fa" {
            let output_path_decompress = output_dir.join("populations.fasta");
            decompress_file(&output_path, &output_path_decompress)?;
        }

        // TBD phylogeny

        Ok(())
    }

    pub fn load(dataset_dir: &Path) -> Result<(), Report> {
        let reference_path = dataset_dir.join("reference.fasta");
        println!("{reference_path:?}");
        Ok(())
    }
}

pub async fn download_file(url: &str, output_path: &PathBuf) -> Result<(), Report> {
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

pub fn decompress_file(input: &PathBuf, output: &PathBuf) -> Result<(), Report> {
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
