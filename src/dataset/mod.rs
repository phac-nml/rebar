use crate::phylogeny::{Phylogeny, PhylogenyExportFormat, PhylogenyImportFormat};
use crate::sequence::{Sequence, Substitution};
use crate::utils;
use crate::ToYaml;
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use color_eyre::section::Section;
use csv;
use itertools::Itertools;
use log::{debug, info};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::default::Default;
use std::fmt;
use std::fs::{create_dir_all, write, File};
use std::path::Path;
use std::str::FromStr;

// ----------------------------------------------------------------------------
// Constants
// ----------------------------------------------------------------------------

pub const SARSCOV2_POPULATIONS_URL: &str = "https://raw.githubusercontent.com/corneliusroemer/pango-sequences/main/data/pango-consensus-sequences_genome-nuc.fasta.zst";
pub const SARSCOV2_REFERENCE_URL: &str = "https://raw.githubusercontent.com/nextstrain/ncov/master/data/references_sequences.fasta";
pub const SARSCOV2_ALIAS_KEY_URL: &str = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json";
pub const SARSCOV2_LINEAGE_NOTES_URL: &str = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt";

// ----------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Dataset Name

#[derive(Copy, Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
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
    pub tag: String,
    pub name: String,
}

impl ToYaml for Summary {}

// ----------------------------------------------------------------------------
// Population Conflict Summary

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ConflictSummary {
    pub support: Vec<Substitution>,
    pub conflict_ref: Vec<Substitution>,
    pub conflict_alt: Vec<Substitution>,
    pub total: isize,
}

impl ConflictSummary {
    pub fn new() -> Self {
        ConflictSummary {
            support: Vec::new(),
            conflict_ref: Vec::new(),
            conflict_alt: Vec::new(),
            total: 0,
        }
    }
}

impl Default for ConflictSummary {
    fn default() -> Self {
        Self::new()
    }
}

// ----------------------------------------------------------------------------
// Dataset Search Result

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SearchResult {
    pub sequence_id: String,
    pub consensus_population: String,
    pub top_populations: Vec<String>,
    pub diagnostic: Vec<String>,
    pub substitutions: Vec<Substitution>,
    pub support: BTreeMap<String, Vec<Substitution>>,
    pub private: Vec<Substitution>,
    pub conflict_ref: BTreeMap<String, Vec<Substitution>>,
    pub conflict_alt: BTreeMap<String, Vec<Substitution>>,
    pub total: BTreeMap<String, isize>,
    pub recombinant: Option<String>,
}

impl ToYaml for SearchResult {
    fn to_yaml(&self) -> String {
        // Order the population lists from 'best' to 'worst'

        // total
        let mut total_order = self.total.iter().collect::<Vec<(_, _)>>();
        total_order.sort_by(|a, b| b.1.cmp(a.1));

        // support
        let mut support_order: Vec<String> = Vec::new();
        for (pop, _count) in &total_order {
            let subs = &self.support[*pop];
            let count = subs.len();
            support_order.push(format!(
                "{}:\n    - count: {}\n    - substitutions: {}",
                pop,
                count,
                subs.iter().join(", ")
            ));
        }

        // conflict_ref
        let mut conflict_ref_order: Vec<String> = Vec::new();
        for (pop, _count) in &total_order {
            let subs = &self.conflict_ref[*pop];
            let count = subs.len();
            conflict_ref_order.push(format!(
                "{}:\n    - count: {}\n    - substitutions: {}",
                pop,
                count,
                subs.iter().join(", ")
            ));
        }

        // conflict_alt
        let mut conflict_alt_order: Vec<String> = Vec::new();
        for (pop, _count) in &total_order {
            let subs = &self.conflict_alt[*pop];
            let count = subs.len();
            conflict_alt_order.push(format!(
                "{}:\n    - count: {}\n    - substitutions: {}",
                pop,
                count,
                subs.iter().join(", ")
            ));
        }

        // Pretty string formatting for yaml
        let total_order = total_order
            .iter()
            .map(|(pop, count)| format!("{}:\n    - count: {}", &pop, &count))
            .collect::<Vec<_>>();

        // private count and list
        let private_order = format!(
            "  - count: {}\n    - substitutions: {}",
            self.private.len(),
            self.private.iter().join(", ")
        );

        format!(
            "sequence_id: {}
consensus_population: {}
top_populations:\n  - {}
diagnostic:\n  - {}
recombinant: {}
substitutions: {}
total:\n  {}
support:\n  {}
conflict_ref:\n  {}
conflict_alt:\n  {}
private:\n  {}",
            self.sequence_id,
            self.consensus_population,
            self.top_populations.join("\n  - "),
            self.diagnostic.join("\n  - "),
            self.recombinant.clone().unwrap_or("None".to_string()),
            self.substitutions.iter().join(", "),
            total_order.join("\n  "),
            support_order.join("\n  "),
            conflict_ref_order.join("\n  "),
            conflict_alt_order.join("\n  "),
            private_order,
        )
    }
}

impl SearchResult {
    pub fn new() -> Self {
        SearchResult {
            sequence_id: String::new(),
            consensus_population: String::new(),
            top_populations: Vec::new(),
            diagnostic: Vec::new(),
            support: BTreeMap::new(),
            private: Vec::new(),
            conflict_ref: BTreeMap::new(),
            conflict_alt: BTreeMap::new(),
            substitutions: Vec::new(),
            total: BTreeMap::new(),
            recombinant: None,
        }
    }
}

impl Default for SearchResult {
    fn default() -> Self {
        Self::new()
    }
}

// ----------------------------------------------------------------------------
// Serializers

#[derive(serde::Deserialize, Debug)]
struct DiagnosticMutationsRow<'a> {
    mutation: &'a str,
    population: &'a str,
    _include_descendants: &'a str,
}

// ----------------------------------------------------------------------------
// Dataset

#[derive(Debug, Deserialize, Serialize)]
pub struct Dataset {
    pub name: Name,
    pub tag: Tag,
    pub reference: Sequence,
    pub populations: BTreeMap<String, Sequence>,
    pub mutations: BTreeMap<Substitution, Vec<String>>,
    pub diagnostic: BTreeMap<Substitution, String>,
    pub phylogeny: Phylogeny,
}

impl fmt::Display for Dataset {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "name: {}, tag: {}", self.name, self.tag)
    }
}

impl Default for Dataset {
    fn default() -> Self {
        Self::new()
    }
}

impl ToYaml for Dataset {}

impl Dataset {
    pub fn new() -> Self {
        Dataset {
            name: Name::Unknown,
            tag: Tag::Unknown,
            reference: Sequence::new(),
            populations: BTreeMap::new(),
            mutations: BTreeMap::new(),
            diagnostic: BTreeMap::new(),
            phylogeny: Phylogeny::new(),
        }
    }

    /// Summarize population conflicts relative to the query sequence.
    pub fn conflict_summary(
        &self,
        population: &String,
        sequence: &Sequence,
        coordinates: &[usize],
    ) -> Result<ConflictSummary, Report> {
        let mut conflict_summary = ConflictSummary::new();

        // get all the substitutions found in this population
        let pop_subs = self.populations[population]
            .substitutions
            .iter()
            .filter(|sub| {
                coordinates.contains(&sub.coord)
                    && !sequence.missing.contains(&sub.coord)
                    && !sequence.deletions.contains(&sub.to_deletion())
            })
            .collect_vec();

        // support: sub in query that is also in pop
        conflict_summary.support = sequence
            .substitutions
            .iter()
            .filter(|sub| coordinates.contains(&sub.coord) && pop_subs.contains(sub))
            .cloned()
            .collect_vec();

        // conflict_alt: sub in query that is not in pop
        conflict_summary.conflict_alt = sequence
            .substitutions
            .iter()
            .filter(|sub| coordinates.contains(&sub.coord) && !pop_subs.contains(sub))
            .cloned()
            .collect_vec();

        // conflict_ref: sub in pop that is not in query sample
        conflict_summary.conflict_ref = pop_subs
            .into_iter()
            .filter(|sub| {
                coordinates.contains(&sub.coord) && !sequence.substitutions.contains(sub)
            })
            .cloned()
            .collect_vec();

        // total: support - conflict_ref
        conflict_summary.total = conflict_summary.support.len() as isize
            - conflict_summary.conflict_ref.len() as isize;

        Ok(conflict_summary)
    }

    /// Search dataset for a population parsimony match to the sequence.
    pub fn search(
        &self,
        sequence: &Sequence,
        populations: Option<&Vec<String>>,
        coordinates: Option<&Vec<usize>>,
    ) -> Result<SearchResult, Report> {
        // initialize an empty result, this will be the final product of this function
        let mut search_result = SearchResult::new();
        search_result.sequence_id = sequence.id.clone();

        // check if we are restricting the population search
        //   otherwise use all populations in the dataset
        let populations_default = self.populations.keys().cloned().collect_vec();
        let populations = match populations {
            Some(pops) => pops,
            None => &populations_default,
        };

        // Check if we are restricting the coordinate search
        // otherwise use all coordinates in genome_length
        let coordinates_default = (1..sequence.genome_length).collect_vec();
        let coordinates = match coordinates {
            Some(coords) => coords,
            None => &coordinates_default,
        };

        // --------------------------------------------------------------------
        // Support

        // check which population have a sub matching the query sequence
        let mut population_matches = Vec::new();

        for sub in &sequence.substitutions {
            // Check if this sub is part of the fn param coordinates
            if !coordinates.contains(&sub.coord) {
                continue;
            }

            // Check if this is a known sub in the database
            if self.mutations.contains_key(sub) {
                // get all populations that have sub
                let matches = self.mutations[sub]
                    .iter()
                    .filter(|pop| {
                        populations.contains(pop) && !population_matches.contains(pop)
                    })
                    .collect_vec();

                // store the matching subs by population
                population_matches.extend(matches);
            } else {
                search_result.private.push(sub.to_owned());
            }
        }

        if population_matches.is_empty() {
            debug!("No mutations matched a population in the dataset.");
            return Ok(search_result);
        }

        // --------------------------------------------------------------------
        // Conflict

        // check which populations have extra subs/lacking subs
        for population in population_matches {
            let conflict_summary =
                self.conflict_summary(population, sequence, coordinates)?;

            // update search results with support and conflicts found
            search_result
                .support
                .insert(population.to_owned(), conflict_summary.support);
            search_result
                .conflict_ref
                .insert(population.to_owned(), conflict_summary.conflict_ref);
            search_result
                .conflict_alt
                .insert(population.to_owned(), conflict_summary.conflict_alt);
            search_result
                .total
                .insert(population.to_owned(), conflict_summary.total);
        }

        // --------------------------------------------------------------------
        // Consensus Population

        // which population(s) has the highest total?
        let max_total = search_result
            .total
            .iter()
            .max_by(|a, b| a.1.cmp(b.1))
            .map(|(_pop, count)| *count)
            .unwrap_or_else(|| {
                panic!("No populations in the summary total for: {}", sequence.id)
            });

        search_result.top_populations = search_result
            .total
            .iter()
            .filter(|(_pop, count)| *count >= &max_total)
            .map(|(pop, _count)| pop)
            .cloned()
            .collect::<Vec<_>>();

        // Undecided if this filter is a good idea
        // But it helps cut down on verbosity and data stored
        search_result.total = search_result
            .total
            .into_iter()
            .filter(|(pop, _count)| search_result.top_populations.contains(pop))
            .collect::<BTreeMap<_, _>>();

        search_result.support = search_result
            .support
            .into_iter()
            .filter(|(pop, _count)| search_result.top_populations.contains(pop))
            .collect::<BTreeMap<_, _>>();

        search_result.conflict_ref = search_result
            .conflict_ref
            .into_iter()
            .filter(|(pop, _subs)| search_result.top_populations.contains(pop))
            .collect::<BTreeMap<_, _>>();

        search_result.conflict_alt = search_result
            .conflict_alt
            .into_iter()
            .filter(|(pop, _count)| search_result.top_populations.contains(pop))
            .collect::<BTreeMap<_, _>>();

        // --------------------------------------------------------------------
        // Outlier Removal

        // check if any diagnostic mutations are present
        // and whether these populations are in the top_populations (?)
        // TBD whether there are cases where the true population
        // does not appear in the top_populations

        search_result.diagnostic = self
            .diagnostic
            .iter()
            .filter(|(sub, pop)| {
                sequence.substitutions.contains(sub)
                    && search_result.top_populations.contains(pop)
            })
            .map(|(_sub, pop)| pop.to_owned())
            .unique()
            .collect::<Vec<_>>();

        // --------------------------------------------------------------------
        // Consensus Population

        // Without a phylogeny, just use first pop in list
        if self.phylogeny.is_empty() {
            // if we found populations with diagnostic mutations, prioritize those
            if !search_result.diagnostic.is_empty() {
                search_result.consensus_population = search_result.diagnostic[0].clone();
            }
            // otherwise use top_populations list
            else {
                search_result.consensus_population =
                    search_result.top_populations[0].clone();
            }
        }
        // Otherwise, summarize top populations by common ancestor
        else {
            // if we found populations with diagnostic mutations, prioritize those
            if !search_result.diagnostic.is_empty() {
                search_result.consensus_population = self
                    .phylogeny
                    .get_common_ancestor(&search_result.diagnostic)?;
            }
            // otherwise use top_populations list
            else {
                search_result.consensus_population = self
                    .phylogeny
                    .get_common_ancestor(&search_result.top_populations)?;
            }
        }

        // if the common_ancestor was not in the populations list, add it
        // this is repeated code from before, maybe make it a function
        if !search_result
            .top_populations
            .contains(&search_result.consensus_population)
        {
            let population = &search_result.consensus_population;
            let conflict_summary =
                self.conflict_summary(population, sequence, coordinates)?;

            // update search results with support and conflicts found
            search_result
                .support
                .insert(population.to_owned(), conflict_summary.support);
            search_result
                .conflict_ref
                .insert(population.to_owned(), conflict_summary.conflict_ref);
            search_result
                .conflict_alt
                .insert(population.to_owned(), conflict_summary.conflict_alt);
            search_result
                .total
                .insert(population.to_owned(), conflict_summary.total);
        }

        // Check if the consensus population is a known recombinant or descendant of one
        for recombinant in self.phylogeny.recombinants.iter() {
            let recombinant_descendants = self.phylogeny.get_descendants(recombinant)?;
            if recombinant_descendants.contains(&search_result.consensus_population) {
                search_result.recombinant = Some(recombinant.to_owned());
            }
        }

        // set consensus population subs
        search_result.substitutions = self.populations
            [&search_result.consensus_population]
            .substitutions
            .iter()
            .filter(|sub| {
                !sequence.missing.contains(&sub.coord)
                    && !sequence.deletions.contains(&sub.to_deletion())
            })
            .map(|sub| sub.to_owned())
            .collect::<Vec<_>>();

        debug!(
            "{}",
            search_result
                .to_yaml()
                .replace('\n', format!("\n{}", " ".repeat(40)).as_str())
        );
        Ok(search_result)
    }
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Dataset Download

/// Download a remote dataset
pub async fn download(name: &Name, tag: &Tag, output_dir: &Path) -> Result<(), Report> {
    if !output_dir.exists() {
        create_dir_all(output_dir)?;
        info!("Creating output directory: {:?}", output_dir);
    }

    // --------------------------------------------------------------------
    // Download Reference

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
    let reference_path = output_dir.join("reference.fasta");
    info!("Downloading reference: {} to {:?}", url, reference_path);
    utils::download_file(&url, &reference_path, decompress).await?;

    // --------------------------------------------------------------------
    // Download Populations

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
    let populations_path = output_dir.join("populations.fasta");
    info!("Downloading populations: {} to {:?}", url, populations_path);
    utils::download_file(&url, &populations_path, decompress).await?;

    // --------------------------------------------------------------------
    // Create Phylogeny

    let phylogeny_path = output_dir.join("phylogeny.dot");
    info!("Creating phylogeny: {:?}", phylogeny_path);

    let mut phylogeny = Phylogeny::new();
    phylogeny.build_graph(name, output_dir).await?;
    // Export to several formats
    phylogeny.export(output_dir, PhylogenyExportFormat::Dot)?;
    phylogeny.export(output_dir, PhylogenyExportFormat::Json)?;

    // --------------------------------------------------------------------
    // Create Summary

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

// ----------------------------------------------------------------------------
// Dataset Load

/// Load a local dataset
pub fn load(dataset_dir: &Path, mask: usize) -> Result<Dataset, Report> {
    // ------------------------------------------------------------------------
    // Load Reference (Required)

    let reference_path = dataset_dir.join("reference.fasta");
    info!("Loading reference: {:?}", reference_path);
    // read in reference from fasta
    let reference_reader =
        fasta::Reader::from_file(reference_path).expect("Unable to load reference");
    // parse just the first record (ie. next() )
    let reference = reference_reader.records().next().unwrap().unwrap();
    // convert to a bio Sequence object
    let reference = Sequence::from_record(reference, None, mask)?;

    // ------------------------------------------------------------------------
    // Load Populations and Parse Mutations (Required)

    let populations_path = dataset_dir.join("populations.fasta");
    info!("Loading populations: {:?}", populations_path);
    // read in reference from fasta
    let populations_reader =
        fasta::Reader::from_file(populations_path).expect("Unable to load populations");
    // store populations and mutations as maps
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

    // ------------------------------------------------------------------------
    // Load Summary (optional)

    let summary_path = dataset_dir.join("summary.yaml");
    let mut tag = Tag::Unknown;
    let mut name = Name::Unknown;
    if summary_path.exists() {
        info!("Loading summary: {:?}", summary_path);
        let reader = File::open(summary_path)?;
        let summary: Summary = serde_yaml::from_reader(&reader)?;
        name = Name::from_str(&summary.name)?;
        tag = std::str::FromStr::from_str(&summary.tag)?;
    }

    // ------------------------------------------------------------------------
    // Load Phylogeny (Optional)
    let phylogeny_path = dataset_dir.join("phylogeny.json");
    let mut phylogeny = Phylogeny::new();
    if phylogeny_path.exists() {
        info!("Loading phylogeny: {:?}", phylogeny_path);
        phylogeny = Phylogeny::import(dataset_dir, PhylogenyImportFormat::Json)?;
    }

    // --------------------------------------------------------------------
    // Load/Create Diagnostic Mutations

    let diagnostic_path = dataset_dir.join("diagnostic_mutations.tsv");
    let mut diagnostic: BTreeMap<Substitution, String> = BTreeMap::new();

    if diagnostic_path.exists() {
        info!("Loading diagnostic mutations.");
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(diagnostic_path)?;
        let _headers = reader.headers()?;
        for result in reader.records() {
            let record = result?;
            let row: DiagnosticMutationsRow = record.deserialize(None)?;
            let mutation = Substitution::from_str(row.mutation)?;
            let population = row.population.to_string();
            diagnostic.insert(mutation.to_owned(), population.clone());
        }
    } else {
        info!("Calculating diagnostic mutations.");
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(diagnostic_path)?;
        let headers = vec!["mutation", "population", "include_descendants"];
        writer.write_record(&headers)?;

        for (mutation, populations) in &mutations {
            let parent = populations[0].to_owned();
            let mut is_diagnostic = false;
            let mut include_descendants = false;

            // If we don't have a phylogeny, a diagnostic mutation is found in only 1 population
            if phylogeny.is_empty() {
                if populations.len() == 1 {
                    is_diagnostic = true;
                }
            }
            // Otherwise with the phylogeny, a diagnostic mutation is monophyletic of the first pop
            else {
                let common_ancestor = phylogeny.get_common_ancestor(populations)?;
                if common_ancestor == parent {
                    is_diagnostic = true;

                    if populations.len() > 1 {
                        include_descendants = true;
                    }
                }
            }

            if is_diagnostic {
                // save to the hashmap
                diagnostic.insert(mutation.to_owned(), parent.clone());
                // write to file as cache
                let row = vec![
                    mutation.to_string(),
                    parent,
                    include_descendants.to_string(),
                ];
                writer.write_record(&row)?;
            }
        }

        writer.flush()?;
    }

    // assembles pieces into full dataset
    let dataset = Dataset {
        name,
        tag,
        reference,
        populations,
        mutations,
        diagnostic,
        phylogeny,
    };

    Ok(dataset)
}
