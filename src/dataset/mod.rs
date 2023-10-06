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

    /// Summarize population conflicts relative to the query sequence.
    pub fn conflict_summary(
        &self,
        population: &String,
        sequence: &Sequence,
        coordinates: Option<&Vec<usize>>,
    ) -> Result<ConflictSummary, Report> {
        let mut conflict_summary = ConflictSummary::new();

        if !self.populations.contains_key(population) {
            return Err(eyre!(
                "Dataset does not contain a sequence for population {population}"
            ));
        }

        // get all the substitutions found in this population
        let mut pop_subs = self.populations[population]
            .substitutions
            .iter()
            .filter(|sub| {
                !sequence.missing.contains(&sub.coord)
                    && !sequence.deletions.contains(&sub.to_deletion())
            })
            .collect_vec();
        if let Some(coordinates) = coordinates {
            pop_subs = pop_subs
                .into_iter()
                .filter(|sub| coordinates.contains(&sub.coord))
                .collect_vec();
        }

        // support: sub in query that is also in pop
        conflict_summary.support = sequence
            .substitutions
            .iter()
            .filter(|sub| pop_subs.contains(sub))
            .cloned()
            .collect_vec();

        // conflict_alt: sub in query that is not in pop
        conflict_summary.conflict_alt = sequence
            .substitutions
            .iter()
            .filter(|sub| !pop_subs.contains(sub))
            .cloned()
            .collect_vec();

        // conflict_ref: sub in pop that is not in query sample
        conflict_summary.conflict_ref = pop_subs
            .into_iter()
            .filter(|sub| !sequence.substitutions.contains(sub))
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
        let mut search_result = SearchResult::new(sequence);

        // check if we are restricting the population search
        //   otherwise use all populations in the dataset
        let populations = if let Some(populations) = populations {
            populations.clone()
        } else {
            self.populations.keys().cloned().collect_vec()
        };

        // --------------------------------------------------------------------
        // Support

        // check which population have a sub matching the query sequence
        let mut population_matches = Vec::new();

        for sub in &sequence.substitutions {
            // Check if this sub is part of the fn param coordinates
            if let Some(coordinates) = coordinates {
                if !coordinates.contains(&sub.coord) {
                    continue;
                }
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
            return Err(eyre!("No mutations matched a population in the dataset."));
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

        debug!("Search Result:\n{}", search_result.pretty_print());
        Ok(search_result)
    }
}

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
    pub total: BTreeMap<String, isize>,
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
            total: BTreeMap::new(),
            recombinant: None,
        }
    }

    fn pretty_print(&self) -> String {
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

        formatdoc!(
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
