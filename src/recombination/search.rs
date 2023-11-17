use crate::cli::run;
use crate::dataset::{Dataset, SearchResult};
use crate::recombination::{detect_recombination, Hypothesis, Recombination};
use crate::sequence::Sequence;
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;
use log::debug;
use std::collections::BTreeMap;
use strum::IntoEnumIterator;

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

/// Search for primary and secondary recombination parents.
///
/// Uses a recursion_limit for safety. It is not intended to
/// run this wrapper function more than once recursively.
pub fn all_parents<'seq>(
    sequence: &'seq Sequence,
    dataset: &Dataset,
    best_match: &SearchResult,
    populations: &[&String],
    args: &run::Args,
) -> Result<Recombination<'seq>, Report> {
    // copy args, we don't want to modify the original global parameters
    let mut args = args.clone();
    // copy search populations, we will refine these based on edge cases and
    // different hypotheses concerning the mechanism of recombination
    let mut populations = populations.to_vec();

    // ------------------------------------------------------------------------
    // Edge Case
    // ------------------------------------------------------------------------

    // // Edge Cases: Manually specified in the organism's dataset.
    let mut edge_case = false;
    let edge_case_search = dataset
        .edge_cases
        .iter()
        .find(|e| e.population.as_ref() == best_match.recombinant.as_ref());

    if let Some(edge_case_args) = edge_case_search {
        debug!("Applying edge case parameters: {edge_case_args:?}");
        edge_case = true;
        args = args.apply_edge_case(edge_case_args)?;

        if let Some(parents) = &edge_case_args.parents {
            populations.retain(|pop| parents.contains(pop))
        }
        if let Some(knockout) = &edge_case_args.knockout {
            populations.retain(|pop| !knockout.contains(pop))
        }
    }

    // ----------------------------------------------------------------------------
    // Hypothesis Testing
    // ----------------------------------------------------------------------------

    // Store the results of our hypothesis testing
    let mut hypotheses: BTreeMap<Hypothesis, (isize, Option<Recombination>)> =
        BTreeMap::new();

    // iterate through the potential hypotheses
    for hypothesis in Hypothesis::iter() {
        debug!("Testing Hypothesis: {hypothesis:?}");

        // ----------------------------------------------------------------------------
        // Hypothesis: Non-Recombinant
        // The sequence is simply it's best match (consensus population)
        // + conflict_alt mutations and + conflict_ref reversions
        if hypothesis == Hypothesis::NonRecombinant {
            if best_match.recombinant.is_none() {
                let parsimony_score =
                    best_match.score.get(&best_match.consensus_population).unwrap();
                hypotheses.insert(Hypothesis::NonRecombinant, (*parsimony_score, None));
            }
            continue;
        }

        let mut hyp_populations: Vec<&String> = populations.clone();

        // ----------------------------------------------------------------------------
        // Hypothesis: Designated Recombinant.
        // The best match (consensus) population is a known recombinant (or descendant of)
        // Search for parents based on the known list of parents.

        if hypothesis == Hypothesis::DesignatedRecombinant {
            if let Some(recombinant) = &best_match.recombinant {
                let designated_parents = dataset.phylogeny.get_parents(recombinant)?;
                debug!("Designated Parents: {designated_parents:?}");
                hyp_populations.retain(|pop| designated_parents.contains(pop));
            }
        }

        // ----------------------------------------------------------------------------
        // Hypothesis: Novel Recombinant, Allow Recursion.
        // One or more parent(s) are recombinants or recombinant descendants
        // This is a default search, no filter edits needed
        if hypothesis == Hypothesis::NovelRecursiveRecombinant {}

        // ----------------------------------------------------------------------------
        // Hypothesis: Novel Recombinant, Disallow Recursion.
        // No parent(s) are recombinants or recombinant descendants
        if hypothesis == Hypothesis::NovelNonRecursiveRecombinant {
            // skip this hypothesis if
            if hypotheses.contains_key(&Hypothesis::NovelNonRecursiveRecombinant) {
                continue;
            }
            hyp_populations
                .retain(|pop| dataset.phylogeny.non_recombinants_all.contains(pop));
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
            if hyp_populations.contains(&&best_match.consensus_population) {
                Ok(best_match.clone())
            } else {
                dataset.search(sequence, Some(&hyp_populations), None)
            };

        // exclude is inverse of
        //let hyp_exclude = dataset.populations.keys().filter(|pop| !hyp_populations.contains(collect_vec()
        let mut hyp_args = args.clone();
        hyp_args.parents = Some(hyp_populations.into_iter().cloned().collect_vec());

        // Check if primary search found anything
        if let Ok(primary_parent) = primary_search {
            debug!("Primary Parent Search was successful.");
            debug!("Secondary Parent(s) Search.");
            let secondary_search = secondary_parents(
                sequence,
                dataset,
                best_match,
                &[primary_parent],
                &hyp_args,
            );

            if let Ok(recombination) = secondary_search {
                debug!("Secondary Parent(s) Search was successful.");
                let parsimony_score: isize = recombination.score.values().sum();

                // adjust the hypothesis, in case it wasn't actually recursive
                let hypothesis = if hypothesis == Hypothesis::NovelRecursiveRecombinant {
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
                        Hypothesis::NovelRecursiveRecombinant
                    } else {
                        Hypothesis::NovelNonRecursiveRecombinant
                    }
                } else {
                    hypothesis
                };

                hypotheses.insert(hypothesis, (parsimony_score, Some(recombination)));
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
                .map(|(hyp, (score, _))| format!("{hyp:?}: {score}"))
                .join(", ")
        );
    }

    // ----------------------------------------------------------------------------
    // Best Hypothesis (highest parsimony score = support - conflict)

    if hypotheses.is_empty() {
        return Err(eyre!("No evidence for any recombination hypotheses."));
    }

    let max_score =
        hypotheses.iter().map(|(_hyp, (score, _recombination))| score).max().unwrap();

    // might want to consider the best_hypothesis as the one with the least conflict
    // rather than best support?
    // Ex. XBB* recursive recombinants, the original recombination (BJ.1 x CJ.1) has
    // super high support, but lots of conflict
    let best_hypothesis = hypotheses
        .iter()
        .filter_map(
            |(hyp, (score, _recombination))| {
                if score == max_score {
                    Some(hyp)
                } else {
                    None
                }
            },
        )
        .next()
        .unwrap();

    if best_hypothesis == &Hypothesis::NonRecombinant {
        return Err(eyre!("Best hypothesis is Non-Recombinant."));
    } else {
        debug!("best_hypothesis: {best_hypothesis:?}");
    }
    // 0: score, 1: recombination
    let mut recombination = hypotheses.get(best_hypothesis).unwrap().1.clone().unwrap();
    recombination.edge_case = edge_case;
    recombination.hypothesis = Some(best_hypothesis.to_owned());

    Ok(recombination)
}

// Search for the secondary recombination parent(s).
pub fn secondary_parents<'seq>(
    sequence: &'seq Sequence,
    dataset: &Dataset,
    best_match: &SearchResult,
    parents: &[SearchResult],
    args: &run::Args,
) -> Result<Recombination<'seq>, Report> {
    // Initialize our 'Recombination' result, that we will modify and update
    // as we iterate through potential parents
    let mut recombination = Recombination::new(sequence);

    // All the parents we know of when the function is starting
    // Will probably just be the primary parent (1)?
    // todo!() test with more than 2 parents
    let mut parents = parents.to_vec();
    let mut num_parents = parents.len();

    // The inclusion parents may have been set by the args
    // Otherwise, use all populations in dataset
    let mut include_populations = if let Some(parents) = &args.parents {
        parents.iter().collect_vec()
    } else {
        dataset.populations.keys().collect_vec()
    };

    // Keep track of how many parent search iterations we've done
    let mut num_iter = 0;

    loop {
        // --------------------------------------------------------------------
        // Loop Break Check: Simple
        // --------------------------------------------------------------------

        // Check if we can break out of the loop based on simple checks like the
        // max number of iterations max number of parents achieved.

        if num_parents >= args.max_parents {
            // Maxing out the number of parents is a SUCCESS
            debug!("Maximum parents reached ({num_parents}).");
            return Ok(recombination);
        }
        if num_iter >= args.max_iter {
            debug!("Maximum iterations reached ({num_iter}).");

            // Finding minimum parents is a SUCCESS
            if num_parents >= args.min_parents {
                return Ok(recombination);
            }
            // otherwise FAILURE
            else {
                return Err(eyre!(
                    "Number of parents ({num_parents}) is less than the minimum ({}).",
                    args.min_parents
                ));
            }
        }

        num_iter += 1;
        debug!("Parent #{}: Iteration {num_iter}", num_parents + 1);

        // --------------------------------------------------------------------
        // Conflict Checks
        // --------------------------------------------------------------------

        // Check for ALT or REF bases present in the sequence that are not resolved
        // by a recombination parent found so far.

        // identify all substitutions found in all parents so far
        let parent_substitutions = parents
            .iter()
            .flat_map(|parent| &parent.substitutions)
            .unique()
            .collect_vec();

        // --------------------------------------------------------------------
        // Conflict ALT

        let conflict_alt = sequence
            .substitutions
            .iter()
            .filter(|sub| !parent_substitutions.contains(sub))
            .collect_vec();

        debug!("conflict_alt: {}", conflict_alt.iter().join(", "));

        // Loop Break Check
        if conflict_alt.len() < args.min_subs {
            debug!("Sufficient conflict_alt resolution reached, stopping parent search.");
            // Finding minimum parents is a SUCCESS
            if num_parents >= args.min_parents {
                return Ok(recombination);
            }
            // Otherwise FAILURE
            else {
                return Err(eyre!(
                    "Number of parents ({num_parents}) is less than the minimum ({}).",
                    args.min_parents
                ));
            }
        }

        // --------------------------------------------------------------------
        // Conflict REF

        // Conflict ref is slightly more complicated, because we don't store
        // REF bases in the dataset. Instead we assume lack of ALT bases means
        // the REF base is present. This is problematic, as it doesn't yet
        // account for indels, missing data, multiallelic sites, etc.

        let conflict_ref = parents
            .iter()
            // get conflict REF between this parent and the sequence
            .flat_map(|p| &p.conflict_ref[&p.consensus_population])
            .unique()
            // search for parents that have the REF base (no ALT)
            .filter(|sub| {
                let parents_with_ref = parents
                    .iter()
                    .filter(|p| !p.substitutions.contains(sub))
                    .collect_vec();
                // No parents have REF base = unresolved
                parents_with_ref.is_empty()
            })
            .collect_vec();
        debug!("conflict_ref: {}", conflict_ref.iter().join(", "));

        // --------------------------------------------------------------------
        // Conflict Combine

        // Combine the conflict REF and conflict ALT coordinates
        // We will focus our parent search on populations that match
        // the query sequence at these coordinates.
        let mut coordinates = conflict_alt
            .iter()
            .map(|sub| sub.coord)
            .chain(conflict_ref.iter().map(|sub| sub.coord))
            .collect_vec();
        coordinates.sort();
        debug!("coordinates: {coordinates:?}");

        if coordinates.is_empty() {
            return Err(eyre!("No coordinates left to search."));
        }

        // --------------------------------------------------------------------
        // EXCLUDE POPULATIONS

        // Remove currently known parents from the search list
        let current_parents =
            parents.iter().map(|result| &result.consensus_population).collect_vec();
        include_populations.retain(|pop| !current_parents.contains(pop));

        // Exclude populations that have substitutions at ALL of the conflict_ref
        if !conflict_ref.is_empty() {
            let conflict_ref_populations = dataset
                .populations
                .iter()
                .filter(|(pop, _seq)| include_populations.contains(pop))
                .filter_map(|(pop, seq)| {
                    let count = seq
                        .substitutions
                        .iter()
                        .filter(|sub| conflict_ref.contains(sub))
                        .collect_vec()
                        .len();
                    (count == conflict_ref.len()).then_some(pop)
                })
                .collect_vec();
            include_populations.retain(|pop| !conflict_ref_populations.contains(pop));
        }
        // --------------------------------------------------------------------
        // INCLUDE POPULATIONS

        // prioritize populations that have min_subs conflict_alt
        // ie. help resolve the conflict_alt

        let conflict_alt_populations = dataset
            .populations
            .iter()
            .filter(|(pop, _seq)| include_populations.contains(pop))
            .filter_map(|(pop, seq)| {
                let count = seq
                    .substitutions
                    .iter()
                    .filter(|sub| conflict_alt.contains(sub))
                    .collect_vec()
                    .len();
                (count >= args.min_subs).then_some(pop)
            })
            .collect_vec();
        include_populations.retain(|pop| conflict_alt_populations.contains(pop));

        // trunclate list for display
        let display_populations = if include_populations.len() <= 10 {
            include_populations.iter().join(", ")
        } else {
            format!("{} ...", include_populations[0..10].iter().join(", "),)
        };
        debug!("Prioritizing conflict_alt resolution: {display_populations:?}");

        // If every single population in the dataset has been excluded, exit here
        // This might because we supplied args.parents that were not actually
        // good hypotheses.
        if include_populations.is_empty() {
            return Err(eyre!("No populations left to search."));
        }

        // --------------------------------------------------------------------
        // Search Dataset #1 and #2 (Coordinate Range vs Precise)
        // --------------------------------------------------------------------

        //let search_modes = vec!["range", "precise"];
        let search_modes = vec!["precise"];
        // In what cases do we want to search the full coordinate range first vs
        // after the specific coordinates?

        for mode in search_modes {
            let search_coords = if mode == "precise" {
                debug!("Searching based on precise coordinates: {coordinates:?}");
                coordinates.clone()
            } else {
                let coord_min = coordinates.iter().min().unwrap();
                let coord_max = coordinates.iter().max().unwrap();
                debug!("Searching based on coordinate range: {coord_min} - {coord_max}");
                (coord_min.to_owned()..coord_max.to_owned()).collect_vec()
            };

            //debug!("dataset.search");

            let parent_candidate = dataset.search(
                sequence,
                Some(&include_populations),
                Some(&search_coords),
            );

            //debug!("detect_recombination");

            // if the search found parents, check for recombination
            if let Ok(parent_candidate) = parent_candidate {
                // remove this parent from future searches
                include_populations
                    .retain(|pop| **pop != parent_candidate.consensus_population);

                // check for recombination
                let detect_result = detect_recombination(
                    sequence,
                    dataset,
                    best_match,
                    &parents,
                    Some(&parent_candidate),
                    &dataset.reference,
                    args,
                );

                // if successful, add this parent to the list and update recombination
                // break out of the search mode loop
                if let Ok(detect_result) = detect_result {
                    num_parents += 1;
                    // reset the iter counter
                    num_iter = 0;
                    parents.push(parent_candidate);
                    recombination = detect_result;
                    break;
                }
            }
        }
    }
}
