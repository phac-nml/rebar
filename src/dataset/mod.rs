pub mod constants;
pub mod name;
pub mod phylogeny;
pub mod summary;
pub mod tag;
pub mod utils;

use crate::dataset::constants::*;
use crate::dataset::name::Name;
use crate::dataset::phylogeny::Phylogeny;
use crate::dataset::summary::Summary;
use crate::dataset::tag::Tag;
use crate::query::match_summary::MatchSummary;
use crate::sequence::{Sequence, Substitution};
use crate::traits::ToYaml;
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, WrapErr};
use color_eyre::section::Section;
use log::info;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::default::Default;
use std::fs::{create_dir_all, write, File};
use std::path::Path;
use std::str::FromStr;

#[derive(Debug, Deserialize, Serialize)]
pub struct Dataset {
    pub name: Name,
    pub tag: Tag,
    pub reference: Sequence,
    pub populations: BTreeMap<String, Sequence>,
    pub mutations: BTreeMap<Substitution, Vec<String>>,
    pub phylogeny: Phylogeny,
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
            mutations: BTreeMap::new(),
            phylogeny: Phylogeny::new(),
        }
    }
}

impl ToYaml for Dataset {}

impl Dataset {
    /// Download a remote dataset
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
        utils::download_file(&url, &output_path, decompress).await?;

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
        utils::download_file(&url, &output_path, decompress).await?;

        // TBD phylogeny

        // --------------------------------------------------------------------
        // Dataset-Specific
        // --------------------------------------------------------------------

        if name == Name::SarsCov2 {
            println!("HERE");
        }
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

    /// Load a local dataset
    pub fn load(dataset_dir: &Path, mask: usize) -> Result<Dataset, Report> {
        // Load the reference (required)
        let reference_path = dataset_dir.join("reference.fasta");
        info!("Loading reference: {:?}", reference_path);
        let reference_reader =
            fasta::Reader::from_file(reference_path).expect("Unable to load reference");
        let reference = reference_reader.records().next().unwrap().unwrap();
        let reference = Sequence::from_record(reference, None, mask)?;

        // Load the populations (required)
        // Also parse mutations
        let populations_path = dataset_dir.join("populations.fasta");
        info!("Loading populations: {:?}", populations_path);
        let populations_reader = fasta::Reader::from_file(populations_path)
            .expect("Unable to load populations");
        let mut populations = BTreeMap::new();
        let mut mutations: BTreeMap<Substitution, Vec<String>> = BTreeMap::new();
        for result in populations_reader.records() {
            let record = result?;
            let sequence = Sequence::from_record(record, Some(&reference), mask)?;
            populations.insert(sequence.id.clone(), sequence.clone());

            for sub in sequence.substitutions {
                mutations
                    .entry(sub)
                    .or_insert(Vec::new())
                    .push(sequence.id.clone());
            }
        }

        // Load the summary (optional)
        let summary_path = dataset_dir.join("summary.yaml");
        let mut tag = Tag::Unknown;
        let mut name = Name::Unknown;
        if summary_path.exists() {
            info!("Loading summary: {:?}", summary_path);
            let reader = File::open(summary_path)?;
            let summary: Summary = serde_yaml::from_reader(&reader)?;
            name = Name::from_str(&summary.name)?;
            tag = Tag::from_str(&summary.tag)?;
        }

        // Load the phylogeny (optional)
        let mut phylogeny = Phylogeny::new();
        let phylogeny_path = dataset_dir.join("phylogeny.dot");
        // Check if it already exists
        if phylogeny_path.exists() {
            info!("Loading phylogeny: {:?}", phylogeny_path);
        }
        // Otherwise we create it
        else {
            info!("Preparing dataset phylogeny: {}", &phylogeny_path.display());
            //phylogeny.build_graph(&name, &tag, &dataset_dir).expect("Failed to build phylogeny.");
        }

        // Finally assemble into dataset collection
        let dataset = Dataset {
            name,
            tag,
            reference,
            populations,
            mutations,
            phylogeny,
        };

        Ok(dataset)
    }

    pub fn find_best_match(
        &self,
        sequence: &Sequence,
        exclude_populations: Option<Vec<String>>,
    ) -> Result<MatchSummary, Report> {
        let mut match_summary = MatchSummary::new();

        // Check if we are excluding certains pop
        let mut exclude_pops = Vec::new();
        if let Some(exclude_populations) = exclude_populations {
            exclude_pops = exclude_populations;
        } else {
            exclude_pops = Vec::new();
        }
        let exclude_populations = exclude_pops;

        // support: check which population have a matching sub
        for sub in &sequence.substitutions {
            if self.mutations.contains_key(sub) {
                let mut matches: Vec<String> = self.mutations[sub]
                    .clone()
                    .iter()
                    .filter(|pop| !exclude_populations.contains(pop))
                    .map(|pop| pop.to_owned())
                    .collect::<Vec<_>>();

                for population in matches {
                    *match_summary
                        .support
                        .entry(population.to_owned())
                        .or_insert(0) += 1;
                }
            } else {
                match_summary.private.push(sub.to_owned());
            }
        }

        if match_summary.support.is_empty() {
            return Ok(match_summary);
        }

        // conflict: check which populations have extra subs/lacking subs
        for population in match_summary.support.keys() {
            let pop_sub = &self.populations[population].substitutions;

            // conflict_ref: sub in pop that is not in query
            let pop_conflict_ref: usize = pop_sub
                .iter()
                .filter(|sub| !sequence.substitutions.contains(sub))
                .collect::<Vec<_>>()
                .len();
            // conflict_alt: sub in query that is not in pop
            let pop_conflict_alt: usize = sequence
                .substitutions
                .iter()
                .filter(|sub| !pop_sub.contains(sub))
                .collect::<Vec<_>>()
                .len();

            let pop_total =
                match_summary.support[population] as isize - pop_conflict_ref as isize;

            match_summary
                .conflict_ref
                .insert(population.to_owned(), pop_conflict_ref);
            match_summary
                .conflict_alt
                .insert(population.to_owned(), pop_conflict_alt);
            match_summary.total.insert(population.to_owned(), pop_total);
        }

        // which population(s) has the highest total?
        let binding = match_summary.total.clone();
        let max_total = binding
            .iter()
            .max_by(|a, b| a.1.cmp(&b.1))
            .map(|(_pop, count)| count)
            .unwrap_or_else(|| {
                panic!("No populations in the summary total for: {}", sequence.id)
            });

        let binding = match_summary.total.clone();
        match_summary.top_populations = binding
            .into_iter()
            .filter(|(_pop, count)| count >= &max_total)
            .map(|(pop, _count)| pop)
            .collect::<Vec<_>>();

        // Undecided if this filter is a good idea
        // But it helps cut down on verbosity and data stored
        match_summary.total = match_summary
            .total
            .into_iter()
            .filter(|(pop, _count)| match_summary.top_populations.contains(pop))
            .collect::<BTreeMap<_, _>>();

        match_summary.support = match_summary
            .support
            .into_iter()
            .filter(|(pop, _count)| match_summary.top_populations.contains(pop))
            .collect::<BTreeMap<_, _>>();

        match_summary.conflict_ref = match_summary
            .conflict_ref
            .into_iter()
            .filter(|(pop, _count)| match_summary.top_populations.contains(pop))
            .collect::<BTreeMap<_, _>>();

        match_summary.conflict_alt = match_summary
            .conflict_alt
            .into_iter()
            .filter(|(pop, _count)| match_summary.top_populations.contains(pop))
            .collect::<BTreeMap<_, _>>();

        // TBD: Remove outliers from top graph
        // Based on diagonstic mutations, phylo distance,

        // TBD: Use graph to get consensus pop
        match_summary.consensus_population = match_summary.top_populations[0].clone();

        Ok(match_summary)
    }
}
