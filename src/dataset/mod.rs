use crate::sequence::Sequence;
use crate::traits::ToYaml;
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, WrapErr};
use color_eyre::section::Section;
use log::info;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::default::Default;
use std::fs::{create_dir_all, remove_file, write, File};
use std::io::Read;
use std::path::{Path, PathBuf};
use tempfile::TempDir;
use zstd::stream::read::Decoder;

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
    Unknown,
}

impl std::fmt::Display for Name {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Name::RsvA => write!(f, "rsv-a"),
            Name::RsvB => write!(f, "rsv-b"),
            Name::SarsCov2 => write!(f, "sars-cov-2"),
            Name::Unknown => write!(f, "unknown"),
        }
    }
}

// Just implement nightly first, we'll figure out latest and archive later
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
pub enum Tag {
    Nightly,
    Latest,
    Archive(String),
    Unknown,
}

impl std::fmt::Display for Tag {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Tag::Latest => write!(f, "latest"),
            Tag::Nightly => write!(f, "nightly"),
            Tag::Archive(tag) => write!(f, "{tag}"),
            Tag::Unknown => write!(f, "unknown"),
        }
    }
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Summary {
    name: String,
    tag: String,
}

impl ToYaml for Summary {}

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

impl Default for Dataset {
    fn default() -> Self {
        Dataset {
            name: Name::Unknown,
            tag: Tag::Unknown,
            reference: Sequence::new(),
            populations: BTreeMap::new(),
        }
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

        let tag = match tag.as_str() {
            "latest" => Tag::Latest,
            "nightly" => Tag::Nightly,
            _ => Tag::Archive(tag.to_string()),
        };

        create_dir_all(output_dir)?;
        info!("Creating output directory: {:?}", output_dir);

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
        let mut decompress = false;
        if ext != "fasta" && ext != "fa" {
            decompress = true;
        }
        let output_path = output_dir.join("reference.fasta");
        info!("Downloading reference: {} to {:?}", url, output_path);
        download_file(&url, &output_path, decompress).await?;

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
        let mut decompress = false;
        if ext != "fasta" && ext != "fa" {
            decompress = true;
        }
        let output_path = output_dir.join("populations.fasta");
        info!("Downloading populations: {} to {:?}", url, output_path);
        download_file(&url, &output_path, decompress).await?;

        // TBD phylogeny

        // --------------------------------------------------------------------
        // Create Summary
        // --------------------------------------------------------------------
        let output_path = output_dir.join("summary.yaml");
        info!("Creating info summary: {:?}", output_path);

        let summary = Summary {
            name: name.to_string(),
            tag: tag.to_string(),
        }
        .to_yaml();

        write(&output_path, summary)
            .wrap_err(format!("Unable to write summary: {:?}", output_path))?;

        Ok(())
    }

    pub fn load(dataset_dir: &Path) -> Result<(), Report> {
        // Load the reference (required)
        let reference_path = dataset_dir.join("reference.fasta");
        info!("Loading reference: {:?}", reference_path);
        let reference_reader =
            fasta::Reader::from_file(reference_path).expect("Unable to load reference");
        let reference = reference_reader.records().next().unwrap().unwrap();
        let reference = Sequence::from_record(reference);

        // Load the populations (required)
        let populations_path = dataset_dir.join("populations.fasta");
        info!("Loading populations: {:?}", populations_path);
        let _populations_reader = fasta::Reader::from_file(populations_path)
            .expect("Unable to load populations");
        // for record in populations_reader.records() {
        //     println!("{record:?}");
        // }

        // Load the summary (optional)
        let summary_path = dataset_dir.join("summary.yaml");
        if summary_path.exists() {
            info!("Loading summary: {:?}", summary_path);
        }

        // Load the phylogeny (optional)
        let phylogeny_path = dataset_dir.join("phylogeny.dot");
        if phylogeny_path.exists() {
            info!("Loading phylogeny: {:?}", phylogeny_path);
        }

        // Debug
        let populations = BTreeMap::new();
        let name = Name::Unknown;
        let tag = Tag::Unknown;
        let _dataset = Dataset {
            name,
            tag,
            reference,
            populations,
        };

        Ok(())
    }
}

pub async fn download_file(
    url: &str,
    output_path: &PathBuf,
    decompress: bool,
) -> Result<(), Report> {
    let ext = Path::new(&url).extension().unwrap().to_str().unwrap();

    let response = reqwest::get(url).await?;
    if response.status() != 200 {
        return Err(eyre!(
            "Unable to download file: {url}\nStatus code {}.",
            response.status()
        ));
    }

    if decompress {
        // Write bytes to a tmp file
        let tmp_dir = TempDir::new()?;
        let tmp_path = PathBuf::from(tmp_dir.path()).join(format!("tmpfile.{ext}"));
        let content = response.bytes().await?;
        write(&tmp_path, content)
            .wrap_err(format!("Unable to write file: {:?}", tmp_path))?;
        decompress_file(&tmp_path, output_path, true)?;
    } else {
        let content = response.text().await?;
        write(output_path, content)
            .wrap_err(format!("Unable to write file: {:?}", output_path))?;
    }

    Ok(())
}

pub fn decompress_file(
    input: &PathBuf,
    output: &PathBuf,
    inplace: bool,
) -> Result<(), Report> {
    let ext = input.extension().unwrap();

    match ext.to_str().unwrap() {
        "zst" => {
            let reader = File::open(input)?;
            let mut decoder = Decoder::new(reader)?;
            let mut buffer = String::new();
            decoder.read_to_string(&mut buffer)?;
            write(output, buffer)
                .wrap_err(format!("Unable to write file: {:?}", output))?;

            if inplace {
                remove_file(input)?;
            }
        }
        _ => return Err(eyre!("Decompression for .{ext:?} is not implemented yet.")),
    };

    Ok(())
}
