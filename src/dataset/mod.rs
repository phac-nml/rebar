pub mod attributes;
pub mod download;
pub mod list;
pub mod load;
pub mod sarscov2;
pub mod toy1;

use crate::cli::run;
use crate::phylogeny::Phylogeny;
use crate::sequence::{parsimony, Sequence, Substitution};
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use indoc::formatdoc;
use itertools::Itertools;
use log::debug;
use noodles::fasta;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::default::Default;
use std::fmt;
use std::fs::File;
use std::io::Write;
use std::path::Path;

// ----------------------------------------------------------------------------
// Dataset

#[derive(Debug, Deserialize, Serialize)]
pub struct Dataset {
    pub name: attributes::Name,
    pub tag: attributes::Tag,
    pub reference: Sequence,
    pub populations: BTreeMap<String, Sequence>,
    pub mutations: BTreeMap<Substitution, Vec<String>>,
    pub phylogeny: Phylogeny,
    pub edge_cases: Vec<run::Args>,
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

impl Dataset {
    pub fn new() -> Self {
        Dataset {
            name: attributes::Name::Custom,
            tag: attributes::Tag::Custom,
            reference: Sequence::new(),
            populations: BTreeMap::new(),
            mutations: BTreeMap::new(),
            phylogeny: Phylogeny::new(),
            edge_cases: Vec::new(),
        }
    }

    pub fn create_consensus(
        &self,
        name: &str,
        populations: &[&str],
    ) -> Result<Sequence, Report> {
        // collect individual population sequences
        let sequences = populations
            .iter()
            .filter_map(|pop| {
                (self.populations.contains_key(*pop)).then_some(&self.populations[*pop])
            })
            .collect_vec();

        // construct consensus
        let consensus = (0..self.reference.genome_length)
            .map(|coord| {
                let bases = sequences.iter().map(|s| s.seq[coord]).unique().collect_vec();
                if bases.len() == 1 {
                    bases[0]
                } else {
                    'N'
                }
            })
            .join("");

        // create new fasta record
        let definition = fasta::record::Definition::new(name, None);
        let sequence = fasta::record::Sequence::from(consensus.as_bytes().to_vec());
        let record = fasta::Record::new(definition, sequence);

        // parse and create Sequence record
        // dataset is already masked, no need
        let mask = Vec::new();

        let sequence = Sequence::from_record(record, Some(&self.reference), &mask)?;

        Ok(sequence)
    }

    /// Expand list of populations with wildcarding.
    pub fn expand_populations(
        &self,
        populations: &[String],
    ) -> Result<Vec<String>, Report> {
        // expand '*' to get descendants
        let expanded = populations
            .iter()
            .map(|p| {
                // if population is '*', use all populations in dataset
                if p == "*" {
                    Ok(self.populations.keys().cloned().collect_vec())
                }
                // if population is 'X*', use all recombinants in dataset
                else if p == "X*" {
                    Ok(self.phylogeny.get_recombinants_all()?)
                }
                // if population ends with '*' expand descendants
                else if p.ends_with('*') {
                    let p = p.replace('*', "");
                    self.phylogeny.get_descendants(&p)
                }
                // simple population name, that is in the dataset
                else if self.populations.contains_key(p) {
                    Ok(vec![p.to_string()])
                } else {
                    Err(eyre!("{p} is not present in the dataset."))
                }
            })
            // flatten and handle the `Result` layer
            .collect::<Result<Vec<_>, Report>>()?
            .into_iter()
            // flatten and handle the `Vec` layer
            .flatten()
            .unique()
            .collect_vec();

        Ok(expanded)
    }

    /// Search dataset for a population parsimony match to the sequence.
    pub fn search(
        &self,
        sequence: &Sequence,
        populations: Option<&Vec<&String>>,
        coordinates: Option<&[usize]>,
    ) -> Result<SearchResult, Report> {
        // initialize an empty result, this will be the final product of this function
        let mut result = SearchResult::new(sequence);

        // --------------------------------------------------------------------
        // Candidate Matches

        // Identify preliminary candidate matches, based on populations that have
        // the greatest number of matching substitutions.
        // NOTE: This is a efficiency shortcut, but the true population is not
        // guaranteed to be in this initial candidate pool.

        // optionally filter subs to the requested coordinates
        let search_subs = if let Some(coordinates) = &coordinates {
            sequence
                .substitutions
                .iter()
                .filter(|sub| coordinates.contains(&sub.coord))
                .collect_vec()
        } else {
            sequence.substitutions.iter().collect()
        };

        // Count up all matching population subs ("support")
        let mut max_support = 0;
        let population_support_counts: BTreeMap<&String, usize> = self
            .populations
            .iter()
            .filter(|(pop, _seq)| {
                if let Some(populations) = populations {
                    populations.contains(pop)
                } else {
                    true
                }
            })
            .filter_map(|(pop, seq)| {
                let count = seq
                    .substitutions
                    .iter()
                    .filter(|sub| search_subs.contains(sub))
                    .collect_vec()
                    .len();
                if count >= max_support {
                    max_support = count
                }
                (count > 0).then_some((pop, count))
            })
            .collect();

        // todo!() decide how much wiggle room we want to give in max support
        // if we want to do max_support - 10, we might need to alter pretty_print
        // so that it only displays the first N candidates (ex. 5,10)
        // this will also cause slow downs
        let population_matches = population_support_counts
            .into_iter()
            .filter_map(|(pop, count)| (count == max_support).then_some(pop))
            //.filter(|(_pop, count)| *count >= (max_support - 10))
            .collect_vec();

        if population_matches.is_empty() {
            return Err(eyre!("No mutations matched a population in the dataset."));
        }

        // --------------------------------------------------------------------
        // Conflict

        // check which populations have extra subs/lacking subs
        population_matches.into_iter().for_each(|pop| {
            // calculate the parsimony score, and store results in map by population
            let pop_seq = &self.populations[pop];
            let summary =
                parsimony::Summary::from_sequence(sequence, pop_seq, coordinates)
                    .unwrap_or_else(|_| {
                        panic!("Failed to create summary from sequence {}", &sequence.id)
                    });
            result.support.insert(pop.to_owned(), summary.support);
            result.conflict_ref.insert(pop.to_owned(), summary.conflict_ref);
            result.conflict_alt.insert(pop.to_owned(), summary.conflict_alt);
            result.score.insert(pop.to_owned(), summary.score);
        });

        // --------------------------------------------------------------------
        // Top Populations
        // Tie breaking, prefer matches with the highest score (support - conflict)
        // beyond that, prefer matches with highest support or lowest conflict?
        // Ex. XCU parent #1 could be FL.23 (highest support) or XBC.1 (lowest conflict)

        // which population(s) has the highest score?
        // reminder: it can be negative when extreme recombinant genomic size
        let max_score =
            result.score.values().max().expect("Failed to get max score of result.");

        let max_score_populations = result
            .score
            .iter()
            .filter_map(|(pop, count)| (count == max_score).then_some(pop))
            .collect_vec();

        // break additional ties by max support
        let max_support = result
            .support
            .iter()
            .filter_map(|(pop, subs)| {
                (max_score_populations.contains(&pop)).then_some(subs.len())
            })
            .max()
            .unwrap_or(0);

        result.top_populations = result
            .score
            .iter()
            .zip(result.support.iter())
            .filter_map(|((pop, score), (_, subs))| {
                (score == max_score && subs.len() == max_support).then_some(pop)
            })
            .cloned()
            .collect_vec();

        // --------------------------------------------------------------------
        // Consensus Population
        // summarize top populations by common ancestor
        let consensus_population = if self.phylogeny.is_empty() {
            //result.top_populations.iter().join("|")
            // just take first?
            result.top_populations[0].clone()
        } else {
            self.phylogeny.get_common_ancestor(&result.top_populations)?
        };
        result.consensus_population = consensus_population.clone();

        // if the common_ancestor was not in the populations list, add it
        let consensus_sequence = if !result
            .top_populations
            .contains(&consensus_population)
        {
            let pop = &consensus_population;

            // // Option #1. Actual sequence of the internal MRCA node?
            // let pop_seq = &self.populations[pop];
            // let summary = parsimony::Summary::from_sequence(sequence, pop_seq, coordinates)?;

            // Option #2. Consensus sequence of top populations?
            let top_populations =
                result.top_populations.iter().map(|s| s.as_ref()).collect_vec();
            debug!("Creating {pop} consensus genome from top populations.");
            let pop_seq = self.create_consensus(pop, &top_populations)?;
            let summary =
                parsimony::Summary::from_sequence(sequence, &pop_seq, coordinates)?;

            // Add consensus summary to search result
            result.support.insert(pop.to_owned(), summary.support);
            result.conflict_ref.insert(pop.to_owned(), summary.conflict_ref);
            result.conflict_alt.insert(pop.to_owned(), summary.conflict_alt);
            result.score.insert(pop.to_owned(), summary.score);

            pop_seq
        } else {
            self
                .populations
                .get(&consensus_population)
                .cloned()
                .unwrap_or_else(|| panic!("Consensus population {consensus_population} is not in the dataset populations."))
        };

        // Filter out non-top populations
        // helps cut down on verbosity in debug log and data stored
        // Ex. XE, lots of BA.2 candidates
        result.score.retain(|p, _| {
            result.top_populations.contains(p) || p == &consensus_population
        });
        result.support.retain(|p, _| {
            result.top_populations.contains(p) || p == &consensus_population
        });
        result.conflict_ref.retain(|p, _| {
            result.top_populations.contains(p) || p == &consensus_population
        });
        result.conflict_alt.retain(|p, _| {
            result.top_populations.contains(p) || p == &consensus_population
        });

        // Check if the consensus population is a known recombinant or descendant of one
        result.recombinant = if self.phylogeny.is_empty() {
            None
        } else {
            self.phylogeny.get_recombinant_ancestor(&consensus_population)?
        };

        // --------------------------------------------------------------------
        // Substitutions
        //  --------------------------------------------------------------------

        // set consensus population subs
        result.substitutions = consensus_sequence
            .substitutions
            .into_iter()
            .filter(|sub| {
                !sequence.missing.contains(&sub.coord)
                    && !sequence.deletions.contains(&sub.to_deletion())
            })
            .collect_vec();

        // private subs (conflict_alt and conflict_ref reversed)
        result.private =
            result.conflict_alt.get(&consensus_population).cloned().unwrap_or_default();
        result
            .conflict_ref
            .get(&consensus_population)
            .cloned()
            .unwrap_or_default()
            .iter()
            .for_each(|sub| {
                let mut sub = *sub;
                std::mem::swap(&mut sub.alt, &mut sub.reference);
                result.private.push(sub);
            });
        result.private.sort();

        debug!("Search Result:\n{}", result.pretty_print());
        Ok(result)
    }

    /// If a population name is in the phylogeny but not in the sequences,
    /// find the closest parent that is in the sequences. Might be itself!
    ///
    /// I don't love this function name, need better!
    pub fn get_ancestor_with_sequence(&self, population: &str) -> Result<String, Report> {
        if self.populations.contains_key(population) {
            return Ok(population.to_string());
        }
        // ancestors can have multiple paths to root, because of recombination
        let ancestors = self.phylogeny.get_ancestors(population)?;
        // filter the ancestor paths to just populations we have sequences for
        // prefer the ancestor path that is the longest
        let ancestors_filter = ancestors
            .into_iter()
            .map(|path| {
                path.into_iter()
                    .filter(|p| self.populations.contains_key(p))
                    .collect_vec()
            })
            .max_by(|a, b| a.len().cmp(&b.len()))
            .unwrap_or_default();

        // use the last element in the path (closest parent)
        let ancestor = ancestors_filter.last();
        if let Some(ancestor) = ancestor {
            Ok(ancestor.clone())
        } else {
            Err(eyre!("No ancestor of {population} has sequence data."))
        }
    }
}

// ----------------------------------------------------------------------------
// Dataset Search Result

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct SearchResult {
    pub sequence_id: String,
    pub consensus_population: String,
    pub top_populations: Vec<String>,
    pub substitutions: Vec<Substitution>,
    pub support: BTreeMap<String, Vec<Substitution>>,
    pub private: Vec<Substitution>,
    pub conflict_ref: BTreeMap<String, Vec<Substitution>>,
    pub conflict_alt: BTreeMap<String, Vec<Substitution>>,
    pub score: BTreeMap<String, isize>,
    pub recombinant: Option<String>,
}

impl SearchResult {
    pub fn new(sequence: &Sequence) -> Self {
        SearchResult {
            sequence_id: sequence.id.clone(),
            consensus_population: String::new(),
            top_populations: Vec::new(),
            support: BTreeMap::new(),
            private: Vec::new(),
            conflict_ref: BTreeMap::new(),
            conflict_alt: BTreeMap::new(),
            substitutions: Vec::new(),
            score: BTreeMap::new(),
            recombinant: None,
        }
    }

    pub fn pretty_print(&self) -> String {
        // Order the population lists from 'best' to 'worst'

        let max_display_items = 10;

        // score
        let mut score_order: Vec<(String, isize)> =
            self.score.clone().into_iter().collect();
        score_order.sort_by(|a, b| b.1.cmp(&a.1));

        // put consensus population first, regardless of score
        let consensus_score: (String, isize) = score_order
            .iter()
            .find(|(pop, _score)| *pop == self.consensus_population)
            .cloned()
            .expect("Failed to order consensus populations by score.");

        score_order.retain(|(pop, _score)| *pop != self.consensus_population);
        score_order.insert(0, consensus_score);

        // restrict display items for brevity
        let display_suffix = if score_order.len() > max_display_items {
            score_order = score_order[0..max_display_items].to_vec();
            "\n  ..."
        } else {
            ""
        };

        let mut support_order: Vec<String> = Vec::new();
        let mut conflict_ref_order: Vec<String> = Vec::new();
        let mut conflict_alt_order: Vec<String> = Vec::new();

        score_order.iter().for_each(|(pop, _count)| {
            let subs = &self.support[pop];
            let count = subs.len();
            support_order.push(format!("- {pop} ({count}): {}", subs.iter().join(", ")));

            let subs = &self.conflict_ref[pop];
            let count = subs.len();
            conflict_ref_order
                .push(format!("- {pop} ({count}): {}", subs.iter().join(", ")));

            let subs = &self.conflict_alt[pop];
            let count = subs.len();
            conflict_alt_order
                .push(format!("- {pop} ({count}): {}", subs.iter().join(", ")));
        });

        // Pretty string formatting for yaml
        let score_order = score_order
            .iter()
            .map(|(pop, count)| format!("- {}: {}", &pop, &count))
            .collect::<Vec<_>>();

        formatdoc!(
            "sequence_id: {}
            consensus_population: {}
            top_populations: {}
            recombinant: {}
            substitutions: {}
            score:\n  {}{display_suffix}
            support:\n  {}{display_suffix}
            conflict_ref:\n  {}{display_suffix}
            conflict_alt:\n  {}{display_suffix}
            private: {}",
            self.sequence_id,
            self.consensus_population,
            self.top_populations.join(", "),
            self.recombinant.clone().unwrap_or("None".to_string()),
            self.substitutions.iter().join(", "),
            score_order.join("\n  "),
            support_order.join("\n  "),
            conflict_ref_order.join("\n  "),
            conflict_alt_order.join("\n  "),
            self.private.iter().join(", ")
        )
    }
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

/// Write mapping of mutations to populations, coordinate sorted.
pub fn write_mutations(
    mutations: &BTreeMap<Substitution, Vec<String>>,
    path: &Path,
) -> Result<(), Report> {
    // convert to vector for coordinate sorting
    let mut mutations = mutations.iter().collect_vec();
    mutations.sort_by(|a, b| a.0.coord.cmp(&b.0.coord));

    // convert substitution to string for serde pretty
    let mutations =
        mutations.iter().map(|(sub, pops)| (sub.to_string(), pops)).collect_vec();
    // create output file
    let mut file = File::create(path)
        .wrap_err_with(|| format!("Failed to create file: {path:?}"))?;

    // parse to string
    let output = serde_json::to_string_pretty(&mutations)
        .wrap_err_with(|| "Failed to parse mutations.".to_string())?;

    // write to file
    file.write_all(format!("{}\n", output).as_bytes())
        .wrap_err_with(|| format!("Failed to write file: {path:?}"))?;

    Ok(())
}
