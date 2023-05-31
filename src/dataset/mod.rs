pub mod constants;
pub mod name;
pub mod phylogeny;
pub mod summary;
pub mod tag;
pub mod utils;

use crate::dataset::constants::*;
use crate::dataset::name::Name;
use crate::dataset::phylogeny::{
    Phylogeny, PhylogenyExportFormat, PhylogenyImportFormat,
};
use crate::dataset::summary::Summary;
use crate::dataset::tag::Tag;
use crate::query::match_summary::MatchSummary;
use crate::sequence::{Sequence, Substitution};
use crate::traits::ToYaml;
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, WrapErr};
use color_eyre::section::Section;
use log::{debug, info};
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

        // --------------------------------------------------------------------
        // Phylogeny
        // --------------------------------------------------------------------

        let output_path = output_dir.join("phylogeny.dot");
        info!("Creating phylogeny: {:?}", output_path);

        let mut phylogeny = Phylogeny::new();
        phylogeny.build_graph(&name, output_dir).await?;
        // Export to several formats
        phylogeny.export(output_dir, PhylogenyExportFormat::Dot)?;
        phylogeny.export(output_dir, PhylogenyExportFormat::Json)?;

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
        // (required) load reference
        let reference_path = dataset_dir.join("reference.fasta");
        info!("Loading reference: {:?}", reference_path);
        let reference_reader =
            fasta::Reader::from_file(reference_path).expect("Unable to load reference");
        let reference = reference_reader.records().next().unwrap().unwrap();
        let reference = Sequence::from_record(reference, None, mask)?;

        // (required) load populations and parse mutations
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

        // (optional) load summary
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

        // (optional) load phylogeny
        let phylogeny_path = dataset_dir.join("phylogeny.json");
        let mut phylogeny = Phylogeny::new();
        if phylogeny_path.exists() {
            info!("Loading phylogeny: {:?}", phylogeny_path);
            phylogeny = Phylogeny::import(dataset_dir, PhylogenyImportFormat::Json)?;
        }

        // assembles pieces into full dataset
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
        exclude_populations: Option<&Vec<String>>,
    ) -> Result<MatchSummary, Report> {
        let mut match_summary = MatchSummary::new();

        // Check if we are excluding certain populations
        let binding = Vec::new();
        let exclude_populations: &Vec<String> = exclude_populations.unwrap_or(&binding);

        // support: check which population have a matching sub
        for sub in &sequence.substitutions {
            if self.mutations.contains_key(sub) {
                let matches: Vec<String> = self.mutations[sub]
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
            .max_by(|a, b| a.1.cmp(b.1))
            .map(|(_pop, count)| count)
            .unwrap_or_else(|| {
                panic!("No populations in the summary total for: {}", sequence.id)
            });

        let binding = match_summary.total.clone();
        match_summary.top_populations = binding
            .into_iter()
            .filter(|(_pop, count)| count >= max_total)
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
        // Based on diagnostic mutations, phylo distance

        // Without a phylogeny, just use first
        if self.phylogeny.is_empty() {
            match_summary.consensus_population = match_summary.top_populations[0].clone();
        }
        // Otherwise, summarize top populations by common ancestor
        else {
            match_summary.consensus_population = self
                .phylogeny
                .get_common_ancestor(&match_summary.top_populations)?;
        }

        // Check if the consensus population is a known recombinant or descendant of one
        for recombinant in self.phylogeny.recombinants.iter() {
            let recombinant_descendants = self.phylogeny.get_descendants(recombinant)?;
            if recombinant_descendants.contains(&match_summary.consensus_population) {
                match_summary.recombinant = Some(recombinant.to_owned());
            }
        }

        debug!(
            "{}",
            match_summary
                .to_yaml()
                .replace('\n', format!("\n{}", " ".repeat(40)).as_str())
        );
        Ok(match_summary)
    }

    pub fn find_parents(
        &self,
        sequence: Sequence,
        best_match: &MatchSummary,
        max_parents: usize,
        max_iter: usize,
    ) -> Result<Vec<MatchSummary>, Report> {
        let mut parents = Vec::new();

        if max_parents == 0 {
            return Ok(parents);
        }

        let mut num_parents = 0;
        let mut exclude_populations = Vec::new();

        // --------------------------------------------------------------------
        // Parent 1
        // --------------------------------------------------------------------

        // Option 1. If this is a known recombinant, exclude the recombinant's
        //   descendants from parent 1 search.
        debug!("parent_1");

        let recombinant = best_match.recombinant.clone();
        let parent_1 = if let Some(recombinant) = recombinant {
            let mut descendants = self.phylogeny.get_descendants(&recombinant)?;
            exclude_populations.append(&mut descendants);
            self.find_best_match(&sequence, Some(&exclude_populations))?
        }
        // Option 2. Not a known recombinant, just use best_match/consensus as parent 1
        else {
            best_match.to_owned()
        };

        // get conflict_ref and conflict_alt for subsequent searches
        let parent_1_subs =
            &self.populations[&parent_1.consensus_population].substitutions;

        // conflict_ref: sub in pop that is not in query
        let conflict_ref = parent_1_subs
            .iter()
            .filter(|sub| !sequence.substitutions.contains(sub))
            .collect::<Vec<_>>();

        // if there are no conflict_ref, this parent is considered a perfect match
        // which indicates little to no evidence of recombination.
        if conflict_ref.is_empty() {
            debug!(
                "No recombination detected, parent_1 ({}) is a perfect match.",
                &parent_1.consensus_population
            );
            return Ok(parents);
        }

        // conflict_alt: sub in query that is not in pop
        let _conflict_alt = sequence
            .substitutions
            .iter()
            .filter(|sub| !parent_1_subs.contains(sub))
            .collect::<Vec<_>>();

        // if there are no conflict_ref, this parent is considered a perfect match
        // which indicates little to no evidence of recombination.
        parents.push(parent_1.clone());
        num_parents += 1;

        let mut num_iter = 0;

        while num_parents < max_parents && num_iter < max_iter {
            num_iter += 1;
            debug!("find_parents iter: {num_iter}");

            // exclude descendants of previous parents
            for parent in &parents {
                let descendants = self
                    .phylogeny
                    .get_descendants(&parent.consensus_population)?;
                let mut descendants_to_add = descendants
                    .iter()
                    .filter(|pop| !exclude_populations.contains(pop))
                    .map(|pop| pop.to_owned())
                    .collect::<Vec<_>>();
                exclude_populations.append(&mut descendants_to_add);
            }

            // we want to keep populations in our search that:
            //   - have at least one conflict_alt
            //   - lack all conflict_ref

            // so exclude_populations will contain populations that:
            //   - have no conflict_alt
            //   - have all conflict_ref
            for sub in &conflict_ref {
                if self.mutations.contains_key(sub) {
                    let populations = &self.mutations[sub];
                    let mut populations_to_exclude = populations
                        .iter()
                        .filter(|pop| !exclude_populations.contains(pop))
                        .map(|pop| pop.to_owned())
                        .collect::<Vec<_>>();

                    exclude_populations.append(&mut populations_to_exclude);
                }
            }

            let match_summary =
                self.find_best_match(&sequence, Some(&exclude_populations))?;

            debug!(
                "parent_{}: {}",
                num_parents + 1,
                match_summary.consensus_population
            );
            debug!(
                "{}",
                match_summary
                    .to_yaml()
                    .replace('\n', format!("\n{}", " ".repeat(40)).as_str())
            );
            // check for robust breakpoints based on cli params

            // ... if breakpoints found

            // Update conflict_ref and conflict_alt based on the parents so far

            num_parents += 1;
        }

        Ok(parents)
    }

    pub fn find_breakpoints(
        &self,
        _sequence: &Sequence,
        _match_summary: &MatchSummary,
    ) -> Result<(), Report> {
        Ok(())
    }
}
