use crate::phylogeny::{Phylogeny, PhylogenyExportFormat, PhylogenyImportFormat};
use crate::recombination::detect_recombination;
use crate::sequence::{Sequence, Substitution};
use crate::traits::ToYaml;
use crate::utils;
use csv;
use bio::io::fasta;
use color_eyre::eyre::{eyre, Result, Report, WrapErr};
use color_eyre::section::Section;
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
// Name
// ----------------------------------------------------------------------------

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
// Tag
// ----------------------------------------------------------------------------

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

// Override the FromStr trait from std
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
// Summary
// ----------------------------------------------------------------------------

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Summary {
    pub tag: String,
    pub name: String,
}

impl ToYaml for Summary {}

// ----------------------------------------------------------------------------
// Search Result
// ----------------------------------------------------------------------------

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SearchResult {
    pub consensus_population: String,
    pub top_populations: Vec<String>,
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
            "consensus_population: {}
top_populations:\n  - {}
recombinant: {}
substitutions: {}
total:\n  {}
support:\n  {}
conflict_ref:\n  {}
conflict_alt:\n  {}
private:\n  {}",
            self.consensus_population,
            self.top_populations.join("\n  - "),
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
            consensus_population: String::new(),
            top_populations: Vec::new(),
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
// ----------------------------------------------------------------------------

#[derive(serde::Deserialize, Debug)]
struct DiagnosticMutationsRow<'a> {
    mutation: &'a str,
    population: &'a str,
    _include_descendants: &'a str,
}

// ----------------------------------------------------------------------------
// Dataset
// ----------------------------------------------------------------------------

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

    /// Search dataset for a population parsimony match to the sequence.
    pub fn search(
        &self,
        sequence: &Sequence,
        populations: Option<&Vec<String>>,
        substitutions: Option<&Vec<Substitution>>,
    ) -> Result<SearchResult, Report> {

        // initialize an empty result, this will be the final product of this function
        let mut search_result = SearchResult::new();
    
        // check if we are restricting the population search
        //   otherwise use all populations in the dataset
        let populations_default = self.populations.keys().cloned().collect_vec();
        let populations = match populations {
            Some(pops) => pops,
            None => &populations_default,
        };
    
        // Check if we are restricting the substitution search
        // otherwise use all subs in query sequence
        let substitutions = match substitutions {
            Some(subs) => subs,
            None => &sequence.substitutions,
        };    
    
        // support: check which population have a sub matching the query sequence
        for sub in substitutions {
            if self.mutations.contains_key(sub){

                // get all populations that have sub
                let mut matches = self.mutations[sub]
                    .iter()
                    .filter(|pop| populations.contains(pop))
                    .collect_vec();
    
                // store the matching subs by population
                for population in matches {
                    search_result
                        .support
                        .entry(population.to_owned())
                        .or_insert(Vec::new())
                        .push(*sub);
                }
            } else {
                search_result.private.push(sub.to_owned());
            }
        }
    
        if search_result.support.is_empty() {
            debug!("No mutations matched a population in the dataset.");
            return Ok(search_result);
        }
    
        // conflict: check which populations have extra subs/lacking subs
        for population in search_result.support.keys() {
            // all the subs in the populations barcode
            let pop_subs = self.populations[population]
                .substitutions
                .iter()
                .filter(|sub| {
                    !sequence.missing.contains(&sub.coord)
                        && !sequence.deletions.contains(&sub.to_deletion())
                })
                .cloned()
                .collect_vec();
    
            // conflict_ref: sub in pop that is not in query
            let pop_conflict_ref = pop_subs
                .iter()
                .filter(|sub| !substitutions.contains(sub))
                .cloned()
                .collect_vec();
    
            // conflict_alt: sub in query that is not in pop
            let pop_conflict_alt = substitutions
                .iter()
                .filter(|sub| !pop_subs.contains(sub))
                .cloned()
                .collect_vec();
    
            let pop_total = search_result.support[population].len() as isize
                - pop_conflict_ref.len() as isize;
    
            search_result
                .conflict_ref
                .insert(population.to_owned(), pop_conflict_ref);
            search_result
                .conflict_alt
                .insert(population.to_owned(), pop_conflict_alt);
            search_result.total.insert(population.to_owned(), pop_total.to_owned());
        }
    
        // which population(s) has the highest total?
        let max_total = search_result.total
            .iter()
            .max_by(|a, b| a.1.cmp(b.1))
            .map(|(_pop, count)| *count)
            .unwrap_or_else(|| {
                panic!("No populations in the summary total for: {}", sequence.id)
            });
    
        search_result.top_populations = search_result.total
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
    
        // check if any diagnostic mutations are present, and whether these populations
        //   are in the top_populations list
        let populations_with_diagnostic_mutations = self.mutations
            .iter()
            .filter(|(sub, _pops)| {
                sequence.substitutions.contains(&sub)
            })
            .map(|(_sub, pops)| pops)
            .flatten()
            .unique()
            .filter(|pop| search_result.top_populations.contains(pop))
            .collect::<Vec<_>>();
    
        // for sub in sequence.substitutions {
        //     if self.mutations.contains_key(&sub) {
        //         let mut populations = self.mutations[&sub];
        //         *populations_with_diagnostic_mutations.extend(&mut populations);
        //     }
        // }
        debug!("diagnostic: {populations_with_diagnostic_mutations:?}");
    
        
    
        // TBD: Remove outliers from top graph
        // Based on diagnostic mutations, phylo distance
    
        // Without a phylogeny, just use first
        if self.phylogeny.is_empty() {
            search_result.consensus_population = search_result.top_populations[0].clone();
        }
        // Otherwise, summarize top populations by common ancestor
        else {
            search_result.consensus_population = self
                .phylogeny
                .get_common_ancestor(&search_result.top_populations)?;
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
    phylogeny.build_graph(&name, output_dir).await?;
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
    // Diagnostic Mutations

    let diagnostic_path = dataset_dir.join("diagnostic_mutations.tsv");
    let mut diagnostic: BTreeMap<Substitution, String> = BTreeMap::new();

    if diagnostic_path.exists(){
        info!("Loading diagnostic mutations.");
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(diagnostic_path)?;
        let _headers = reader.headers()?;
        for result in reader.records() {
            let record = result?;
            let row: DiagnosticMutationsRow = record.deserialize(None)?;
            let mutation = Substitution::from_str(row.mutation)?;
            let population =  row.population.to_string();
            diagnostic.insert(mutation.to_owned(), population.clone());
        }
    }
    else {

        info!("Calculating diagnostic mutations.");
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(diagnostic_path)?; 
        let headers = vec!("mutation", "population", "include_descendants");
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
                let common_ancestor = phylogeny.get_common_ancestor(&populations)?;
                if common_ancestor == parent {
                    is_diagnostic = true;

                    if populations.len()  > 1 {
                        include_descendants = true;
                    }
                }
            }

            if is_diagnostic {
                // save to the hashmap
                diagnostic.insert(mutation.to_owned(), parent.clone());
                // write to file as cache
                let row = vec!(mutation.to_string(), parent, include_descendants.to_string());
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

#[allow(clippy::too_many_arguments)]
pub fn find_parents(
    dataset: &Dataset,
    sequence: Sequence,
    best_match: &SearchResult,
    max_parents: usize,
    max_iter: usize,
    min_consecutive: usize,
    min_length: usize,
    min_subs: usize,
) -> Result<Vec<SearchResult>, Report> {
    let mut parents = Vec::new();

    if max_parents == 0 {
        return Ok(parents);
    }

    let mut num_parents = 0;
    // by default, include all dataset populations
    let mut include_populations = dataset.populations.keys().map(|pop| pop.to_owned()).collect::<Vec<_>>();
    // by default, include all sequence substitutions
    let mut include_substitutions = sequence.substitutions.clone();
    // by default, don't exclude any subs or populations
    let mut exclude_populations = Vec::new();
    let mut exclude_substitutions: Vec<Substitution> = Vec::new();

    // --------------------------------------------------------------------
    // Parent 1
    // --------------------------------------------------------------------

    // Option 1. If this is a known recombinant, exclude the recombinant's
    //   descendants from parent 1 search.

    debug!("parent_1");
    let recombinant = best_match.recombinant.clone();
    let parent_search_result = if let Some(recombinant) = recombinant {
        let descendants = dataset.phylogeny.get_descendants(&recombinant)?;
        exclude_populations.extend(descendants);
        include_populations = include_populations
            .iter()
            .filter(|pop| !exclude_populations.contains(&pop))
            .map(|pop| pop.to_owned())
            .collect::<Vec<_>>();
        dataset.search(&sequence, Some(&include_populations), Some(&include_substitutions))?
    }
    // Option 2. Not a known recombinant, just use best_match/consensus as parent 1
    else {
        best_match.to_owned()
    };

    parents.push(parent_search_result);
    num_parents += 1;

   // --------------------------------------------------------------------
    // Parents 2-MAX
    // --------------------------------------------------------------------

    let mut num_iter = 0;

    loop {

        // --------------------------------------------------------------------
        // Loop Break Checks

        if num_parents >= max_parents {
            debug!("Maximum number of parents reached ({max_parents}), stopping parent search.");
            break;
        }
        if num_iter >= max_iter {
            debug!("Maximum number of iterations reached ({num_iter}), stopping parent search.");
            break;
        }

        num_iter += 1;

        // --------------------------------------------------------------------
        // Current Mutation Conflicts

        // first identify all substitutions in all parents
        let mut parent_substitutions = parents
            .iter()
            .map(|p| p.substitutions.clone())
            .flatten()
            .collect_vec();
        parent_substitutions = parent_substitutions
            .iter()
            .unique()
            .cloned()
            .collect_vec();

        // identify conflict_alt that are not resolved by another other parent
        let mut conflict_alt = parents
            .iter()
            .map(|p| p.conflict_alt[&p.consensus_population].clone())
            .flatten()
            .filter(|sub| !parent_substitutions.contains(&sub))
            .collect_vec();
        conflict_alt = conflict_alt
            .iter()
            .unique()
            .cloned()
            .collect_vec();

        // conflict_ref is strange, if any parent DOESN'T HAVE the sub, that means it's resolved
        // not entirely ideal, what about missing and deletions?
        let mut conflict_ref = parents
            .iter()
            .map(|p| p.conflict_ref[&p.consensus_population].clone())
            .flatten()
            .collect_vec();
        conflict_ref = conflict_ref
            .iter()
            .unique()
            .cloned()
            .collect_vec();

        let mut conflict_ref_resolved = Vec::new();
        for sub in &conflict_ref {
            let is_resolved = false;
            for parent in &parents {
                if !parent.substitutions.contains(sub)
                {
                    conflict_ref_resolved.push(sub.to_owned());
                }
            }
        }

        conflict_ref = conflict_ref.into_iter().filter(|sub| !conflict_ref_resolved.contains(sub)).collect_vec();     

        // --------------------------------------------------------------------
        // Loop Break Checks

        if conflict_ref.is_empty() {
            debug!("Sufficient conflict_ref resolution reached, stopping parent search.");
            break;
        }
        if conflict_alt.len() < min_subs {
            debug!("Sufficient conflict_alt resolution reached, stopping parent search.");
            break;
        }                

        debug!("parent_{}: iteration {num_iter}", num_parents + 1);
        debug!("conflict_ref: {}", conflict_ref.iter().join(", ")); 
        debug!("conflict_alt: {}", conflict_alt.iter().join(", "));

        // --------------------------------------------------------------------
        // Current Mutation Support

        let mut parent_support = parents
            .iter()
            .map(|p| p.support[&p.consensus_population].clone())
            .flatten()
            .collect_vec();
        parent_support = parent_support
            .iter()
            .unique()
            .cloned()
            .collect_vec();
        debug!("sequence_resolved: {}", parent_support.iter().join(", "));

        let resolved = sequence.substitutions.iter().filter(|sub| !parent_support.contains(&sub)).cloned().collect_vec();
        debug!("sequence_unresolved: {}", resolved.iter().join(", "));

        // --------------------------------------------------------------------
        // Filters (include/exclude)

        // exclude descendants and ancestors(?) of all previous parents
        for parent in &parents {

            // descendants to exclude
            let descendants = dataset
                .phylogeny
                .get_descendants(&parent.consensus_population)?;
            let descendants_to_exclude = descendants
                .into_iter()
                .filter(|pop| !exclude_populations.contains(pop))
                .collect::<Vec<_>>();
            exclude_populations.extend(descendants_to_exclude);

            // ancestors to exclude
            let ancestors = dataset
                .phylogeny
                .get_ancestors(&parent.consensus_population)?;
            // because of possible recombination, ancestors is a vector of vectors
            // to reflect multiple parents and paths to the root.
            // just flatten them for all our purposes here
            let ancestors_to_add = ancestors
                .into_iter()
                .flatten()
                .filter(|pop| !exclude_populations.contains(pop))
                .collect::<Vec<_>>();
            exclude_populations.extend(ancestors_to_add);            
        }

        // we want to EXCLUDE populations in our search that
        //   - have ALL of the conflict_ref
        let mut conflict_ref_count = BTreeMap::new();
        for sub in &conflict_ref {
            if dataset.mutations.contains_key(sub) {
                let populations = &dataset.mutations[sub];
                for pop in populations {
                    *conflict_ref_count.entry(pop).or_insert(0) += 1
                }
            }            
        }
        let populations_to_exclude = conflict_ref_count
            .into_iter()
            .filter(|(pop, count)| {
                *count == conflict_ref.len() 
                && !exclude_populations.contains(pop)
                && !include_populations.contains(pop) 
            })
            .map(|(pop, _count)| pop.to_owned())
            .collect::<Vec<_>>();        
        exclude_populations.extend(populations_to_exclude);

        // we want to INCLUDE populations in our search that:
        //   - have AT LEAST min_subs conflict_alt

        // count up the number of conflict_alt by population
        let mut conflict_alt_count = BTreeMap::new();
        for sub in &conflict_alt {
            if dataset.mutations.contains_key(sub) {
                let populations = &dataset.mutations[sub];
                for pop in populations {
                    *conflict_alt_count.entry(pop).or_insert(0) += 1
                }
            }            
        }

        // reset the include list to just these populations
        include_populations = conflict_alt_count
            .into_iter()
            .filter(|(pop, count)| {
                *count >= min_subs
                && !exclude_populations.contains(pop)
            })
            .map(|(pop, _count)| pop.to_owned())
            .collect::<Vec<_>>();

        // exclude substitutions that we've already found support for?
        //include_substitutions = include_substitutions.iter().filter(|sub| )

        // // reset the include substitutions list
        // debug!("conflict_ref: {conflict_ref:?}");
        // include_substitutions = conflict_ref.iter().cloned().cloned().collect::<Vec<_>>();
        // include_substitutions.extend(conflict_alt);
        // // include_substitutions = conflict_ref.extend(conflict_alt);

        // DEBUG
        //include_populations = vec!(String::from("BA.2.75"));

        let search_result = dataset.search(
            &sequence,
            Some(&include_populations),
            Some(&include_substitutions),
        )?;
        let recombination = detect_recombination(
            &sequence,
            &parents,
            &search_result,
            min_consecutive,
            min_length,
            min_subs,
        )?;

        // if the recombination search failed, exclude search_result top_populations from next iteration
        if recombination.breakpoints.is_empty() {
            let mut populations_to_exclude = search_result
                .top_populations
                .iter()
                .filter(|pop| !exclude_populations.contains(pop))
                .map(|pop| pop.to_owned())
                .collect::<Vec<_>>();
            exclude_populations.append(&mut populations_to_exclude);
        } else {
            num_parents += 1;
            parents.push(search_result);
        }
    }

    Ok(parents)
}
