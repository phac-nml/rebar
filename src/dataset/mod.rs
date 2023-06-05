use crate::phylogeny::{Phylogeny, PhylogenyExportFormat, PhylogenyImportFormat};
use crate::query::match_summary::MatchSummary;
use crate::recombination::detect_recombination;
use crate::sequence::{Sequence, Substitution};
use crate::traits::ToYaml;
use crate::utils;
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, WrapErr};
use color_eyre::section::Section;
use log::{debug, info};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::default::Default;
use std::fmt;
use std::fs::{create_dir_all, write, File};
use std::path::Path;
use std::str::FromStr;

// ----------------------------------------------------------------------------
// Name
// ----------------------------------------------------------------------------

pub const SARSCOV2_POPULATIONS_URL: &str = "https://raw.githubusercontent.com/corneliusroemer/pango-sequences/main/data/pango-consensus-sequences_genome-nuc.fasta.zst";
pub const SARSCOV2_REFERENCE_URL: &str = "https://raw.githubusercontent.com/nextstrain/ncov/master/data/references_sequences.fasta";
pub const SARSCOV2_ALIAS_KEY_URL: &str = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json";
pub const SARSCOV2_LINEAGE_NOTES_URL: &str = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt";

// ----------------------------------------------------------------------------
// Name
// ----------------------------------------------------------------------------

#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
pub enum Name {
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

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub enum Tag {
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
            // _ => Err(eyre!("Unknown dataset tag: {tag}")).suggestion("Please choose from:")?,
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
// Dataset
// ----------------------------------------------------------------------------

#[derive(Debug, Deserialize, Serialize)]
pub struct Dataset {
    pub name: Name,
    pub tag: Tag,
    pub reference: Sequence,
    pub populations: BTreeMap<String, Sequence>,
    pub mutations: BTreeMap<Substitution, Vec<String>>,
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
            phylogeny: Phylogeny::new(),
        }
    }
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

/// Download a remote dataset
pub async fn download(name: &str, tag: &str, output_dir: &Path) -> Result<(), Report> {
    let name = match name {
        "rsv-a" | "rsv-b" => {
            Err(eyre!("Dataset download of {name} is not implemented yet."))?
        }
        "sars-cov-2" => Name::SarsCov2,
        _ => Err(eyre!("Unknown dataset name: {name}"))
            .suggestion("Please choose from:")?,
    };

    let tag = match tag {
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
    let populations_reader =
        fasta::Reader::from_file(populations_path).expect("Unable to load populations");
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
    dataset: &Dataset,
    sequence: &Sequence,
    exclude_populations: Option<&Vec<String>>,
    include_populations: Option<&Vec<String>>,
) -> Result<MatchSummary, Report> {
    let mut match_summary = MatchSummary::new();

    // Check if we are excluding/including certain populations
    let exclude_binding = Vec::new();
    let exclude_populations: &Vec<String> =
        exclude_populations.unwrap_or(&exclude_binding);

    let include_binding = Vec::new();
    let include_populations: &Vec<String> =
        include_populations.unwrap_or(&include_binding);

    // support: check which population have a matching sub
    for sub in &sequence.substitutions {
        if dataset.mutations.contains_key(sub) {
            let mut matches: Vec<String> = dataset.mutations[sub]
                .clone()
                .iter()
                .filter(|pop| !exclude_populations.contains(pop))
                .map(|pop| pop.to_owned())
                .collect::<Vec<_>>();

            if !include_populations.is_empty() {
                matches = matches
                    .into_iter()
                    .filter(|pop| include_populations.contains(pop))
                    .collect::<Vec<_>>();
            }

            for population in matches {
                match_summary
                    .support
                    .entry(population.to_owned())
                    .or_insert(Vec::new())
                    .push(*sub);
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
        let pop_subs = &dataset.populations[population]
            .substitutions
            .iter()
            .filter(|sub| {
                !sequence.missing.contains(&sub.coord)
                    && !sequence.deletions.contains(&sub.to_deletion())
            })
            .map(|sub| sub.to_owned())
            .collect::<Vec<_>>();

        // conflict_ref: sub in pop that is not in query
        let pop_conflict_ref = pop_subs
            .iter()
            .filter(|sub| !sequence.substitutions.contains(sub))
            .map(|sub| sub.to_owned())
            .collect::<Vec<_>>();

        // conflict_alt: sub in query that is not in pop
        let pop_conflict_alt = sequence
            .substitutions
            .iter()
            .filter(|sub| !pop_subs.contains(sub))
            .map(|sub| sub.to_owned())
            .collect::<Vec<_>>();

        let pop_total = match_summary.support[population].len() as isize
            - pop_conflict_ref.len() as isize;

        match_summary
            .conflict_ref
            .insert(population.to_owned(), pop_conflict_ref.to_owned());
        match_summary
            .conflict_alt
            .insert(population.to_owned(), pop_conflict_alt.to_owned());
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
        .filter(|(pop, _subs)| match_summary.top_populations.contains(pop))
        .collect::<BTreeMap<_, _>>();

    match_summary.conflict_alt = match_summary
        .conflict_alt
        .into_iter()
        .filter(|(pop, _count)| match_summary.top_populations.contains(pop))
        .collect::<BTreeMap<_, _>>();

    // TBD: Remove outliers from top graph
    // Based on diagnostic mutations, phylo distance

    // Without a phylogeny, just use first
    if dataset.phylogeny.is_empty() {
        match_summary.consensus_population = match_summary.top_populations[0].clone();
    }
    // Otherwise, summarize top populations by common ancestor
    else {
        match_summary.consensus_population = dataset
            .phylogeny
            .get_common_ancestor(&match_summary.top_populations)?;
    }

    // Check if the consensus population is a known recombinant or descendant of one
    for recombinant in dataset.phylogeny.recombinants.iter() {
        let recombinant_descendants = dataset.phylogeny.get_descendants(recombinant)?;
        if recombinant_descendants.contains(&match_summary.consensus_population) {
            match_summary.recombinant = Some(recombinant.to_owned());
        }
    }

    // set consensus population subs
    match_summary.substitutions = dataset.populations
        [&match_summary.consensus_population]
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
        match_summary
            .to_yaml()
            .replace('\n', format!("\n{}", " ".repeat(40)).as_str())
    );
    Ok(match_summary)
}

#[allow(clippy::too_many_arguments)]
pub fn find_parents(
    dataset: &Dataset,
    sequence: Sequence,
    best_match: &MatchSummary,
    max_parents: usize,
    max_iter: usize,
    min_consecutive: usize,
    min_length: usize,
    min_subs: usize,
) -> Result<Vec<MatchSummary>, Report> {
    let mut parents = Vec::new();

    if max_parents == 0 {
        return Ok(parents);
    }

    let mut num_parents = 0;
    let mut include_populations = Vec::new();
    let mut exclude_populations = Vec::new();

    // --------------------------------------------------------------------
    // Parent 1
    // --------------------------------------------------------------------

    // Option 1. If this is a known recombinant, exclude the recombinant's
    //   descendants from parent 1 search.

    debug!("parent_1");
    let recombinant = best_match.recombinant.clone();
    let parent_1 = if let Some(recombinant) = recombinant {
        let mut descendants = dataset.phylogeny.get_descendants(&recombinant)?;
        exclude_populations.append(&mut descendants);
        find_best_match(
            dataset,
            &sequence,
            Some(&exclude_populations),
            Some(&include_populations),
        )?
    }
    // Option 2. Not a known recombinant, just use best_match/consensus as parent 1
    else {
        best_match.to_owned()
    };

    // get conflict_ref and conflict_alt for subsequent searches
    let parent_1_subs =
        &dataset.populations[&parent_1.consensus_population].substitutions;

    // conflict_ref: sub in pop that is not in query
    let conflict_ref = parent_1_subs
        .iter()
        .filter(|sub| !sequence.substitutions.contains(sub))
        .collect::<Vec<_>>();

    // if there are no conflict_ref, this parent is considered a perfect match
    // which indicates little to no evidence of recombination.
    if conflict_ref.is_empty() {
        debug!(
            "{} is a perfect match, stopping parent search.",
            &parent_1.consensus_population
        );
        return Ok(parents);
    }

    // conflict_alt: sub in query that is not in pop
    let conflict_alt = sequence
        .substitutions
        .iter()
        .filter(|sub| !parent_1_subs.contains(sub))
        .collect::<Vec<_>>();

    // if there are no conflict_ref, this parent is considered a perfect match
    // which indicates little to no evidence of recombination.
    parents.push(parent_1.clone());
    num_parents += 1;

    let mut num_iter = 0;

    loop {
        if num_parents >= max_parents {
            debug!("Maximum number of parents reached ({max_parents}), stopping parent search.");
            break;
        }
        if num_iter >= max_iter {
            debug!("Maximum number of iterations reached ({num_iter}), stopping parent search.");
            break;
        }

        num_iter += 1;

        // exclude descendants of previous parents
        for parent in &parents {
            let descendants = dataset
                .phylogeny
                .get_descendants(&parent.consensus_population)?;
            let mut descendants_to_add = descendants
                .iter()
                .filter(|pop| !exclude_populations.contains(pop))
                .map(|pop| pop.to_owned())
                .collect::<Vec<_>>();
            exclude_populations.append(&mut descendants_to_add);
        }

        // we want to include populations in our search that:
        //   - have at least one conflict_alt
        //   - are not already in the exclude list

        for sub in &conflict_alt {
            if dataset.mutations.contains_key(sub) {
                let populations_to_add = dataset.mutations[sub]
                    .iter()
                    .filter(|pop| {
                        !exclude_populations.contains(pop)
                            && !include_populations.contains(pop)
                    })
                    .map(|pop| pop.to_owned())
                    .collect::<Vec<_>>();
                include_populations.extend(populations_to_add);
            }
        }

        // we want to exclude populations in our search that
        //    - have all the conflict_ref

        // for sub in &conflict_ref {
        //     if dataset.mutations.contains_key(sub) {
        //         let populations = &dataset.mutations[sub];
        //         for pop in populations {
        //             *conflict_ref_count.entry(pop.clone()).or_insert(0) += 1
        //         }
        //     }
        // }

        debug!("include_populations: {include_populations:?}");

        // // filter out conflict_ref populations
        // for (pop, count) in conflict_ref_count {
        //     if count == conflict_ref.len() && !exclude_populations.contains(&pop) {
        //         exclude_populations.push(pop);
        //     }
        // }
        // // filter out conflict_alt populations
        // for (pop, rcount) in conflict_alt_count {
        //     if alt_count == conflict_alt.len() && alt_count == 0 && !exclude_populations.contains(&pop) {
        //         exclude_populations.push(pop);
        //     }
        // }

        debug!("exclude_populations: {exclude_populations:?}");

        debug!("parent_{}: iteration {num_iter}", num_parents + 1);
        let match_summary = find_best_match(
            dataset,
            &sequence,
            Some(&exclude_populations),
            Some(&include_populations),
        )?;
        let recombination = detect_recombination(
            &sequence,
            &parents,
            &match_summary,
            min_consecutive,
            min_length,
            min_subs,
        )?;

        // if the recombination search failed, exclude match_summary descendants from next iteration
        if recombination.breakpoints.is_empty() {
            let descendants = dataset
                .phylogeny
                .get_descendants(&match_summary.consensus_population)?;
            let mut descendants_to_add = descendants
                .iter()
                .filter(|pop| !exclude_populations.contains(pop))
                .map(|pop| pop.to_owned())
                .collect::<Vec<_>>();
            exclude_populations.append(&mut descendants_to_add);
        } else {
            num_parents += 1;

            // TBD: Check if all conflict_ref are resolved
        }
    }

    Ok(parents)
}
