use crate::cli::run;
use crate::dataset::{Dataset, SearchResult};
use crate::recombination::{detect_recombination, validate, Hypothesis, Recombination};
use crate::sequence::Sequence;
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;
use log::{debug, warn};
use std::collections::BTreeMap;
use strum::IntoEnumIterator;

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

/// Search for primary and secondary recombination parents.
///
/// Uses a recursion_limit for safety. It is not intended to
/// run this wrapper function more than once recursively.
#[allow(clippy::needless_if)]
pub fn all_parents(
    sequence: &Sequence,
    dataset: &Dataset,
    best_match: &mut SearchResult,
    populations: &[String],
    args: &run::Args,
) -> Result<Recombination, Report> {
    // copy args, we don't want to modify the original global parameters
    let mut args = args.clone();
    // copy search populations, we will refine these based on edge cases and
    // different hypotheses concerning the mechanism of recombination
    let mut populations = populations.to_vec();

    let consensus_population = &best_match.consensus_population;

    // ------------------------------------------------------------------------
    // Edge Case
    // ------------------------------------------------------------------------

    // // Edge Cases: Manually specified in the organism's dataset.
    let mut edge_case = false;

    // only apply edge cases if the user didn't request a naive search
    if !args.naive {
        let edge_case_search = dataset
            .edge_cases
            .iter()
            .find(|e| e.population.as_ref() == best_match.recombinant.as_ref());

        if let Some(edge_case_args) = edge_case_search {
            debug!("Applying edge case parameters: {edge_case_args:?}");
            edge_case = true;
            args = args.apply_edge_case(edge_case_args)?;

            if let Some(parents) = &edge_case_args.parents {
                let parents_expand = dataset.expand_populations(parents)?;
                populations.retain(|pop| parents_expand.contains(pop));
            }
            if let Some(knockout) = &edge_case_args.knockout {
                let knockout = dataset.expand_populations(knockout)?;
                populations.retain(|pop| !knockout.contains(pop));
            }
        }
    }

    // ----------------------------------------------------------------------------
    // Hypothesis Testing
    // ----------------------------------------------------------------------------

    // Store the results of our hypothesis testing
    // Hyp: (Recombination, Parents, score, conflict)
    let mut hypotheses: BTreeMap<
        Hypothesis,
        (Option<Recombination>, Vec<SearchResult>, isize, usize),
    > = BTreeMap::new();

    // iterate through the potential hypotheses
    for hypothesis in Hypothesis::iter() {
        debug!("Testing Hypothesis: {hypothesis:?}");

        // ----------------------------------------------------------------------------
        // Hypothesis: Non-Recombinant
        // The sequence is simply it's best match (consensus population)
        // + conflict_alt mutations and + conflict_ref reversions
        if hypothesis == Hypothesis::NonRecombinant {
            if best_match.recombinant.is_none() {
                let score = best_match.score[consensus_population];
                let conflict = best_match.conflict_alt[consensus_population].len()
                    + best_match.conflict_ref[consensus_population].len();
                hypotheses.insert(
                    Hypothesis::NonRecombinant,
                    (None, Vec::new(), score, conflict),
                );
            }
            continue;
        }

        let mut hyp_populations = populations.clone();

        // ----------------------------------------------------------------------------
        // Hypothesis: Designated Recombinant.
        // The best match (consensus) population is a known recombinant (or descendant of)
        // Search for parents based on the known list of parents.

        if hypothesis == Hypothesis::DesignatedRecombinant {
            // skip this hypothesis if user requested naive search
            if args.naive || best_match.recombinant.is_none() {
                continue;
            }
            if let Some(recombinant) = &best_match.recombinant {
                let mut designated_parents =
                    dataset.phylogeny.get_parents(recombinant)?;
                debug!("Designated parents: {designated_parents:?}");

                // we might not have sequence data for all designated parents.
                let designated_parents_filter = designated_parents
                    .iter()
                    .filter_map(|p| dataset.get_ancestor_with_sequence(p).ok())
                    .collect_vec();

                if designated_parents != designated_parents_filter {
                    debug!("Designated parents with sequence data: {designated_parents_filter:?}");
                    // if the filter yielded a different length, one parent is unsearchable
                    if designated_parents.len() != designated_parents_filter.len() {
                        debug!("DesignatedRecombinant hypothesis is impossible to test! At least one parent has no ancestor with sequence data!");
                        continue;
                    }
                    designated_parents = designated_parents_filter;
                }

                // Option #1. no wildcarding, parental strict
                // hyp_populations.retain(|pop| designated_parents.contains(pop));

                // Option #2. add wildcarding for all descendants (exclude novel recombination)
                // ex. XJ is designated as BA.1 and BA.2, specifically, BD.1 (BA.1.17.2.1) and BA.2.65
                // ex  XDF is designated as "XBB*","EG.5.1.3", specifically EG.10.1 and EG.5.1.3

                let designated_wildcards = designated_parents
                    .into_iter()
                    .map(|p| format!("{p}*-r"))
                    .collect_vec();
                let designated_parents_expanded = dataset
                    .expand_populations(&designated_wildcards)?
                    .into_iter()
                    .unique()
                    .collect_vec();

                hyp_populations.retain(|pop| designated_parents_expanded.contains(pop));

                // exclude recombinant itself and descendants from parent search
                //let recombinant_descendants = dataset.phylogeny.get_descendants(recombinant)?;
                //hyp_populations.retain(|pop| !recombinant_descendants.contains(pop));
            }
        }

        // ----------------------------------------------------------------------------
        // Hypothesis: Recombinant, Allow Recursion.
        // One or more parent(s) are recombinants or recombinant descendants
        // This is a default search, no filter edits needed
        if hypothesis == Hypothesis::RecursiveRecombinant {}

        // ----------------------------------------------------------------------------
        // Hypothesis: Recombinant, Allow Recursion, Knockout
        // Don't allow best_match to be primary parent
        // One or more parent(s) are recombinants or recombinant descendants
        // This is a default search, no filter edits needed
        if hypothesis == Hypothesis::KnockoutRecombinant {
            hyp_populations.retain(|pop| *pop != best_match.consensus_population);
        }

        // ----------------------------------------------------------------------------
        // Hypothesis: Non Recursive Recombinant, Disallow Recursion.
        // No parent(s) are recombinants or recombinant descendants
        if hypothesis == Hypothesis::NonRecursiveRecombinant {
            // skip this hypothesis if
            if hypotheses.contains_key(&Hypothesis::NonRecursiveRecombinant) {
                continue;
            }
            hyp_populations
                .retain(|pop| !dataset.phylogeny.recombinants_all.contains(pop));
        }

        // ----------------------------------------------------------------------------
        // Hypothesis Test Time!
        // Search for primary and scondary parents

        debug!("Primary Parent Search.");
        if hyp_populations.is_empty() {
            debug!("No parent populations fit this hypothesis.");
            continue;
        }

        // we only need to rerun the primary parent search if the search populations
        // DOES NOT contain the hyp_populations
        let primary_search: Result<SearchResult, Report> =
            if hyp_populations.contains(&best_match.consensus_population) {
                Ok(best_match.clone())
            } else {
                dataset.search(sequence, Some(&hyp_populations), None)
            };

        // exclude is inverse of
        //let hyp_exclude = dataset.populations.keys().filter(|pop| !hyp_populations.contains(collect_vec()
        let mut hyp_args = args.clone();
        hyp_args.parents = Some(hyp_populations);

        // Check if primary search found anything
        if let Ok(primary_parent) = primary_search {
            debug!("Primary Parent Search was successful.");
            debug!("Secondary Parent(s) Search.");
            let secondary_search =
                secondary_parents(sequence, dataset, &[primary_parent], &hyp_args);

            if let Ok((recombination, parents)) = secondary_search {
                debug!("Secondary Parent(s) Search was successful.");
                let score: isize = recombination.score.values().sum();
                let conflict_alt: usize =
                    recombination.conflict_alt.values().map(|subs| subs.len()).sum();
                let conflict_ref: usize =
                    recombination.conflict_ref.values().map(|subs| subs.len()).sum();
                let conflict = conflict_alt + conflict_ref;

                // adjust the hypothesis, in case it wasn't actually recursive
                let hypothesis = if hypothesis == Hypothesis::RecursiveRecombinant {
                    let mut is_recursive = false;
                    recombination.parents.iter().for_each(|pop| {
                        let recombinant_ancestor = dataset
                            .phylogeny
                            .get_recombinant_ancestor(pop)
                            .unwrap_or(None);
                        if recombinant_ancestor.is_some() {
                            is_recursive = true
                        }
                    });

                    if is_recursive {
                        Hypothesis::RecursiveRecombinant
                    } else {
                        Hypothesis::NonRecursiveRecombinant
                    }
                } else {
                    hypothesis
                };

                hypotheses
                    .insert(hypothesis, (Some(recombination), parents, score, conflict));
            } else {
                debug!("Secondary Parent(s) Search was unsuccessful.");
            }
        } else {
            debug!("Primary Parent Search was unsuccessful.");
        }

        debug!(
            "Hypotheses: {}",
            hypotheses
                .iter()
                .map(|(hyp, (_r, _p, score, conflict))| format!(
                    "{hyp:?}: score={score}, conflict={conflict}"
                ))
                .join(", ")
        );
    }

    // ----------------------------------------------------------------------------
    // Best Hypothesis

    // No evidence at all, return error
    if hypotheses.is_empty() {
        return Err(eyre!("No evidence for any recombination hypotheses."));
    }
    // single hypothesis
    let best_hypothesis = if hypotheses.len() == 1 {
        hypotheses.first_key_value().unwrap().0.clone()
    }
    // evidence for multiple hypotheses, take max score
    else {
        let max_score = hypotheses.iter().map(|(_hyp, (_r, _p, s, _c))| s).max().unwrap();
        let best_hypotheses = hypotheses
            .iter()
            .filter_map(|(hyp, (_r, _p, s, _c))| (s == max_score).then_some(hyp))
            .collect_vec();

        // if hypotheses are tied, prefer them in enum order (first before last)
        // note: this currently means Designated is preferred over non-designated.
        //       prefer this for now, since we can use --naive to disable designated
        let hypothesis_ranks: BTreeMap<Hypothesis, usize> =
            Hypothesis::iter().enumerate().map(|(i, h)| (h, i)).collect();
        let best_hyp_rank = best_hypotheses
            .iter()
            .map(|hyp| {
                hypothesis_ranks
                    .get(hyp)
                    .expect("Hypothesis ranks does not contain hypothesis {hyp:?}")
            })
            .min()
            .unwrap();
        let best_hypothesis = hypothesis_ranks
            .iter()
            .filter_map(|(hyp, r)| (r == best_hyp_rank).then_some(hyp))
            .next()
            .unwrap();

        best_hypothesis.to_owned()
    };

    debug!("best_hypothesis: {best_hypothesis:?}");

    // non-recombinant means that the parent search "failed"
    if best_hypothesis == Hypothesis::NonRecombinant {
        return Err(eyre!("Best hypothesis is Non-Recombinant."));
    }
    let result = hypotheses
        .remove(&best_hypothesis)
        .expect("Hypotheses does not contain the best hypothesis {best_hypothesis:?}");
    let mut recombination = result.0.unwrap();
    let primary_parent = result.1[0].clone();

    // ------------------------------------------------------------------------
    // Recombinant attributes

    recombination.edge_case = edge_case;

    // Decide on novel vs known recombinant at this point
    recombination.recombinant = if let Some(recombinant) = &best_match.recombinant {
        // check if expected parents match observed
        let observed = &recombination.parents;
        let expected = dataset.phylogeny.get_parents(recombinant)?;
        let parents_match = validate::compare_parents(observed, &expected, dataset)?;

        if parents_match {
            Some(recombinant.clone())
        } else {
            *best_match = primary_parent;
            debug!(
                "Changing best match consensus to primary parent:\n{}",
                best_match.pretty_print()
            );
            Some("novel".to_string())
        }
    } else {
        Some("novel".to_string())
    };

    recombination.unique_key = format!(
        "{}_{}_{}",
        &recombination.recombinant.clone().unwrap(),
        &recombination.parents.iter().join("_"),
        &recombination.breakpoints.iter().join("_"),
    );

    recombination.hypothesis = Some(best_hypothesis);

    Ok(recombination)
}

// Search for the secondary recombination parent(s).
pub fn secondary_parents(
    sequence: &Sequence,
    dataset: &Dataset,
    parents: &[SearchResult],
    args: &run::Args,
) -> Result<(Recombination, Vec<SearchResult>), Report> {
    // Initialize our 'Recombination' result, that we will modify and update
    // as we iterate through potential parents
    let mut recombination = Recombination::new();

    // All the parents we know of when the function is starting
    // Will probably just be the primary parent (1)?
    // todo!() test with more than 2 parents
    let mut parents = parents.to_vec();
    let mut num_parents = parents.len();

    // The inclusion parents may have been set by the args
    // Otherwise, use all populations in dataset
    let mut include_populations = if let Some(parents) = &args.parents {
        parents.iter().cloned().collect_vec()
    } else {
        dataset.populations.keys().cloned().collect_vec()
    };

    // Keep track of how many parent search iterations we've done
    let mut num_iter = 0;

    loop {
        // --------------------------------------------------------------------
        // Loop Break Check: Simple
        // --------------------------------------------------------------------

        // Check if we can break out of the loop based on simple checks like the
        // max number of iterations max number of par