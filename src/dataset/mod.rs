pub mod attributes;
pub mod download;
pub mod list;
pub mod load;
//pub mod io;
pub mod sarscov2;

use crate::cli::run;
use crate::phylogeny::Phylogeny;
use crate::sequence::{Sequence, Substitution};
use color_eyre::eyre::{eyre, Report, Result};
use indoc::formatdoc;
use itertools::Itertools;
use log::debug;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::default::Default;
use std::fmt;

// ----------------------------------------------------------------------------
// Dataset

#[derive(Debug, Deserialize, Serialize)]
pub struct Dataset {
    pub name: attributes::Name,
    pub tag: attributes::Tag,
    pub reference: Sequence,
    pub populations: BTreeMap<String, Sequence>,
    pub mutations: BTreeMap<Substitution, Vec<String>>,
    pub diagnostic: BTreeMap<Substitution, String>,
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
            name: attributes::Name::Unknown,
            tag: attributes::Tag::Unknown,
            reference: Sequence::new(),
            populations: BTreeMap::new(),
            mutations: BTreeMap::new(),
            diagnostic: BTreeMap::new(),
            phylogeny: Phylogeny::new(),
            edge_cases: Vec::new(),
        }
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

    /// Summarize population support and conflicts relative to the query sequence.
    pub fn parsimony_summary(
        &self,
        population: &String,
        sequence: &Sequence,
        coordinates: Option<&Vec<usize>>,
    ) -> Result<ParsimonySummary, Report> {
        let mut parsimony_summary = ParsimonySummary::new();

        if !self.populations.contains_key(population) {
            return Err(eyre!(
                "Dataset does not contain a sequence for population {population}"
            ));
        }

        // get all the substitutions found in this population
        // exclude missing and deletion coordinates
        let mut pop_subs = self.populations[population]
            .substitutions
            .iter()
            .filter(|sub| {
                !sequence.missing.contains(&sub.coord)
                    && !sequence.deletions.contains(&sub.to_deletion())
            })
            .collect_vec();

        // get all the substitutions found in the sequence
        let mut seq_subs = sequence.substitutions.clone();

        // optionally filter coordinates
        if let Some(coordinates) = coordinates {
            pop_subs.retain(|sub| coordinates.contains(&sub.coord));
            seq_subs.retain(|sub| coordinates.contains(&sub.coord));
        }

        // support: sub in query seq that is also in pop
        // conflict_alt: sub in query seq that is not in pop
        seq_subs.iter().for_each(|sub| {
            if pop_subs.contains(&sub) {
                parsimony_summary.support.push(*sub);
            } else {
                parsimony_summary.conflict_alt.push(*sub);
            }
        });

        // conflict_ref: sub in pop that is not in query seq
        parsimony_summary.conflict_ref = pop_subs
            .into_iter()
            .filter(|sub| !seq_subs.contains(sub))
            .cloned()
            .collect_vec();

        // score: support - conflict_alt - conflict_ref
        // why did we previously use only conflict_ref and not conflict_alt?
        // individual isize conversion otherwise: "attempt to subtract with overflow"
        parsimony_summary.score = parsimony_summary.support.len() as isize
            - parsimony_summary.conflict_ref.len() as isize
            - parsimony_summary.conflict_alt.len() as isize;

        Ok(parsimony_summary)
    }

    /// Search dataset for a population parsimony match to the sequence.
    pub fn search(
        &self,
        sequence: &Sequence,
        populations: Option<&Vec<String>>,
        coordinates: Option<&Vec<usize>>,
    ) -> Result<SearchResult, Report> {
        // initialize an empty result, this will be the final product of this function
        let mut search_result = SearchResult::new(sequence);

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
                .cloned()
                .collect_vec()
        } else {
            sequence.substitutions.clone()
        };

        // Count up all matching population subs ("support")
        let mut population_support_counts = BTreeMap::new();
        search_subs
            .iter()
            .filter(|sub| self.mutations.contains_key(sub))
            .for_each(|sub| {
                self.mutations[sub]
                    .iter()
                    .filter(|pop| {
                        // optionally filter populations
                        if let Some(populations) = &populations {
                            populations.contains(pop)
                        } else {
                            true
                        }
                    })
                    .for_each(|pop| {
                        population_support_counts
                            .entry(pop)
                            .and_modify(|count| *count += 1)
                            .or_insert(0);
                    })
            });

        // Get maximum support value
        let max_support = population_support_counts
            .values()
            .cloned()
            .max()
            .unwrap_or(0);

        // todo!() decide how much wiggle room we want to give in max support
        // if we want to do max_support - 10, we might need to alter pretty_print
        // so that it only displays the first N candidates (ex. 5,10)
        let population_matches = population_support_counts
            .into_iter()
            //.filter(|(_pop, count)| *count == max_support)
            .filter(|(_pop, count)| *count >= (max_support - 10))
            .map(|(pop, _count)| pop)
            .cloned()
            .collect_vec();

        if population_matches.is_empty() {
            return Err(eyre!("No mutations matched a population in the dataset."));
        }

        // --------------------------------------------------------------------
        // Conflict

        // check which populations have extra subs/lacking subs
        population_matches.iter().for_each(|pop| {
            let summary = self.parsimony_summary(pop, sequence, coordinates).unwrap();
            search_result
                .support
                .insert(pop.to_owned(), summary.support);
            search_result
                .conflict_ref
                .insert(pop.to_owned(), summary.conflict_ref);
            search_result
                .conflict_alt
                .insert(pop.to_owned(), summary.conflict_alt);
            search_result.score.insert(pop.to_owned(), summary.score);
        });

        // --------------------------------------------------------------------
        // Top Populations
        // Tie breaking, prefer matches with the highest score (support - conflict)
        // beyond that, prefer matches with highest support or lowest conflict?
        // Ex. XCU parent #1 could be FL.23 (highest support) or XBC.1 (lowest conflict)

        // which population(s) has the highest score?
        let max_score = search_result.score.values().max().unwrap();

        search_result.top_populations = search_result
            .score
            .iter()
            .filter(|(_pop, count)| *count == max_score)
            .map(|(pop, _count)| pop)
            .cloned()
            .collect_vec();

        // break additional ties by max support
        let max_support = search_result
            .support
            .iter()
            .filter(|(pop, _subs)| search_result.top_populations.contains(pop))
            .map(|(_pop, subs)| subs.len())
            .max()
            .unwrap();

        search_result.top_populations = search_result
            .score
            .iter()
            .zip(search_result.support.iter())
            .filter(|((_, score), (_, subs))| {
                *score == max_score && subs.len() == max_support
            })
            .map(|((pop, _score), _)| pop)
            .cloned()
            .collect_vec();

        // --------------------------------------------------------------------
        // Diagnostic Mutations

        // todo!()
        // check if any diagnostic mutations are present
        // and whether these populations are in the top_populations (?)
        // todo!()
        //  whether there are cases where the true population
        // does not appear in the top_populations

        // search_result.diagnostic = self
        //     .diagnostic
        //     .iter()
        //     .filter(|(sub, pop)| {
        //         sequence.substitutions.contains(sub)
        //             && search_result.top_populations.contains(pop)
        //     })
        //     .map(|(_sub, pop)| pop.to_owned())
        //     .unique()
        //     .collect::<Vec<_>>();

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
            let pop = &search_result.consensus_population;
            let summary = self.parsimony_summary(pop, sequence, coordinates)?;

            // update search results with support and conflicts found
            search_result
                .support
                .insert(pop.to_owned(), summary.support);
            search_result
                .conflict_ref
                .insert(pop.to_owned(), summary.conflict_ref);
            search_result
                .conflict_alt
                .insert(pop.to_owned(), summary.conflict_alt);
            search_result.score.insert(pop.to_owned(), summary.score);
        }

        // Check if the consensus population is a known recombinant or descendant of one
        search_result.recombinant = self
            .phylogeny
            .get_recombinant_ancestor(&search_result.consensus_population)?;

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

        // Filter out non-top populations from the debug log.
        // helps cut down on verbosity and data stored
        // Ex. XE, lots of BA.2 candidates
        search_result.score = search_result
            .score
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

        debug!("Search Result:\n{}", search_result.pretty_print());
        Ok(search_result)
    }
}

// ----------------------------------------------------------------------------
// Population Parsimony Summary

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ParsimonySummary {
    pub support: Vec<Substitution>,
    pub conflict_ref: Vec<Substitution>,
    pub conflict_alt: Vec<Substitution>,
    pub score: isize,
}

impl ParsimonySummary {
    pub fn new() -> Self {
        ParsimonySummary {
            support: Vec::new(),
            conflict_ref: Vec::new(),
            conflict_alt: Vec::new(),
            score: 0,
        }
    }

    pub fn pretty_print(&self) -> String {
        formatdoc!(
            "score:\n  {}
            support:\n  {}
            conflict_ref:\n  {}
            conflict_alt:\n  {}",
            self.score,
            self.support.iter().join(", "),
            self.conflict_ref.iter().join(", "),
            self.conflict_alt.iter().join(", "),
        )
    }
}

impl Default for ParsimonySummary {
    fn default() -> Self {
        Self::new()
    }
}

// ----------------------------------------------------------------------------
// Dataset Search Result

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
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
    pub score: BTreeMap<String, isize>,
    pub recombinant: Option<String>,
}

impl SearchResult {
    pub fn new(sequence: &Sequence) -> Self {
        SearchResult {
            sequence_id: sequence.id.clone(),
            consensus_population: String::new(),
            top_populations: Vec::new(),
            diagnostic: Vec::new(),
            support: BTreeMap::new(),
            private: Vec::new(),
            conflict_ref: BTreeMap::new(),
            conflict_alt: BTreeMap::new(),
            substitutions: Vec::new(),
            score: BTreeMap::new(),
            recombinant: None,
        }
    }

    fn pretty_print(&self) -> String {
        // Order the population lists from 'best' to 'worst'

        // score
        let mut score_order = self.score.iter().collect::<Vec<(_, _)>>();
        score_order.sort_by(|a, b| b.1.cmp(a.1));

        // support
        let mut support_order: Vec<String> = Vec::new();
        for (pop, _count) in &score_order {
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
        for (pop, _count) in &score_order {
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
        for (pop, _count) in &score_order {
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
        let score_order = score_order
            .iter()
            .map(|(pop, count)| format!("{}:\n    - count: {}", &pop, &count))
            .collect::<Vec<_>>();

        // private count and list
        let private_order = format!(
            "  - count: {}\n    - substitutions: {}",
            self.private.len(),
            self.private.iter().join(", ")
        );

        formatdoc!(
            "sequence_id: {}
            consensus_population: {}
            top_populations:\n  - {}
            diagnostic:\n  - {}
            recombinant: {}
            substitutions: {}
            score:\n  {}
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
            score_order.join("\n  "),
            support_order.join("\n  "),
            conflict_ref_order.join("\n  "),
            conflict_alt_order.join("\n  "),
            private_order,
        )
    }
}
