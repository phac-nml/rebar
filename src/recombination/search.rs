use crate::cli::run;
use crate::dataset::{Dataset, SearchResult};
use crate::recombination::{detect_recombination, Recombination};
use crate::sequence::Sequence;
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;
use log::debug;
use std::collections::BTreeMap;

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
    args: &run::Args,
) -> Result<(Vec<SearchResult>, Recombination<'seq>), Report> {
    let mut args = args.clone();
    let mut exclude_populations: Vec<String> = Vec::new();
    let mut include_populations = dataset.populations.keys().cloned().collect_vec();
    let mut edge_case = true;
    let best_match_pop = best_match.consensus_population.to_string();

    // ------------------------------------------------------------------------
    // Edge Case : Low Priority

    let edge_case_args = dataset
        .edge_cases
        .iter()
        .find(|e| e.population.as_ref() == Some(&best_match_pop));

    if let Some(edge_case_args) = edge_case_args {
        debug!("Applying edge case parameters: {edge_case_args:?}");
        edge_case = true;
        args = args.apply_edge_case(edge_case_args)?;

        if let Some(populations) = &edge_case_args.parents {
            include_populations = populations.clone();
        }
        if let Some(populations) = &edge_case_args.knockout {
            exclude_populations = populations.clone();
        }
    }

    // ------------------------------------------------------------------------
    // args.parents : Medium Priority
    // args.parents supplied on the CLI overrides edge case

    if let Some(populations) = &args.parents {
        include_populations = populations.clone();
    }

    // ------------------------------------------------------------------------
    // args.knockout : High Priority
    // args.knockout supplied on the CLI overrides edge cases and args.parents

    if let Some(populations) = &args.knockout {
        exclude_populations = populations.clone();
    }

    include_populations.retain(|p| !exclude_populations.contains(p));
    exclude_populations.retain(|p| !include_populations.contains(p));

    // ----------------------------------------------------------------------------
    // Primary Parent Search

    let mut parents = Vec::new();

    debug!("Primary Parent Search.");
    let parent_primary = dataset.search(sequence, Some(&include_populations), None)?;
    let parent_primary_pop = parent_primary.consensus_population.to_string();
    parents.push(parent_primary);

    // ----------------------------------------------------------------------------
    // Secondary Parents
    debug!("Secondary Parents Search.");

    // Decide whether consensus descendants can be secondary parents
    if !include_populations.contains(&best_match_pop) {
        debug!("Removing best match {best_match_pop}* from search.");
        let descendants = dataset.phylogeny.get_descendants(&best_match_pop)?;
        exclude_populations.extend(descendants);
    }

    // Decide whether parent #1 descendants can be secondary parents
    if !include_populations.contains(&parent_primary_pop) {
        debug!("Removing parent #1 {parent_primary_pop}* from search.");
        let descendants = dataset.phylogeny.get_descendants(&parent_primary_pop)?;
        exclude_populations.extend(descendants);
    }

    let mut result = secondary_parents(
        sequence,
        dataset,
        best_match,
        &parents,
        &exclude_populations,
        &args,
    );

    // If first time failed, try again, this time don't let the consensus
    // population be a parent
    if result.is_err() && best_match_pop == parent_primary_pop {
        debug!(
            "Attempting another search, best match {best_match_pop} cannot be parent."
        );

        // Add the best match and its descendants to the knockout
        let best_match_descendants =
            dataset.phylogeny.get_descendants(&best_match_pop)?;
        args.knockout = if let Some(mut populations) = args.knockout {
            populations.extend(best_match_descendants);
            Some(populations)
        } else {
            Some(best_match_descendants)
        };

        result = all_parents(sequence, dataset, best_match, &args);
    }

    let (parents, mut recombination) = result?;

    recombination.edge_case = edge_case;

    Ok((parents, recombination))
}
// // ----------------------------------------------------------------------------
// // Check known vs novel status

// // If there was more than 1 parent, but the best match is not a known recombinant
// // this is a novel recombinant.
// if parents.len() > 1 && best_match.recombinant.is_none() {
//     best_match.recombinant = Some("novel".to_string());
// }
// // check if observed parents match expected parents
// // ex. XBL (before designated) was XBB ( BA.2.75 and XBB.1.5.57)
// // which conflicts with expected XBB parents BJ.1 (BA.2.10.1) and BM.1.1.1 (BA.2.75.3)
// // todo!() decide on the order, should descendants be from observed or expected?
// else if let Some(recombinant) = &best_match.recombinant {

//     let expected_parents = dataset.phylogeny.get_parents(&recombinant)?;
//     let mut observed_parents_descendants = recombination.parents
//         .iter()
//         .flat_map(|p| dataset.phylogeny.get_descendants(p).unwrap())
//         .unique()
//         .collect_vec();

//     observed_parents_descendants.retain(|p| expected_parents.contains(p));
//     if observed_parents_descendants.len() != expected_parents.len(){
//         best_match.recombinant = Some("novel".to_string());
//     }
// }

// Search for the secondary recombination parent(s).
pub fn secondary_parents<'seq>(
    sequence: &'seq Sequence,
    dataset: &Dataset,
    best_match: &SearchResult,
    parents: &[SearchResult],
    exclude_populations: &[String],
    args: &run::Args,
) -> Result<(Vec<SearchResult>, Recombination<'seq>), Report> {
    // Initialize our 'Recombination' result, that we will modify and update
    // as we iterate through potential parents
    let mut recombination = Recombination::new(sequence);

    // All the parents we know of when the function is starting
    // Will probably just be the primary parent (1)?
    let mut parents = parents.to_vec();
    let mut num_parents = parents.len();

    // Control which populations we will/will not examine as candidate parents
    let mut exclude_populations = exclude_populations.to_vec();
    // The inclusion parents may have been set by the args
    let mut include_populations = if let Some(parents) = &args.parents {
        parents.clone()
    } else {
        Vec::new()
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
            return Ok((parents, recombination));
        }
        if num_iter >= args.max_iter {
            debug!("Maximum iterations reached ({num_iter}).");

            // Finding minimum parents is a SUCCESS
            if num_parents >= args.min_parents {
                return Ok((parents, recombination));
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
                return Ok((parents, recombination));
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
        // the REF base is present. This is problematic, as it doesn't account
        // for indels or missing data yet.

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

        // No loop break checks for conflict_ref.
        // todo!() add examples of why we don't break on conflict_ref resolution.

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

        // --------------------------------------------------------------------
        // EXCLUDE POPULATIONS

        // Add currently known parents to search exclusion
        // Don't exclude their descendants! (ex. XBL is BA.2.75 and XBB.1.5.57)
        for p in &parents {
            if exclude_populations.contains(&p.consensus_population) {
                exclude_populations.push(p.consensus_population.clone());
            }
        }
        // Exclude populations that have substitutions at ALL of the conflict_ref
        let mut conflict_ref_count = BTreeMap::new();
        conflict_ref
            .iter()
            // Make sure this mutation is in the dataset (not private)
            .filter(|sub| dataset.mutations.contains_key(sub))
            // Identify populations that have ALT base
            .for_each(|sub| {
                dataset.mutations[sub]
                    .iter()
                    // Keep track of how many ALT bases each population has
                    .for_each(|p| *conflict_ref_count.entry(p).or_insert(0) += 1);
            });
        for (p, count) in conflict_ref_count.into_iter() {
            if !exclude_populations.contains(p) && count == conflict_ref.len() {
                exclude_populations.push(p.to_string());
            }
        }

        // --------------------------------------------------------------------
        // INCLUDE POPULATIONS

        include_populations.retain(|p| !exclude_populations.contains(p));

        // If the include_populations list is empty at this point
        // ex. No parents supplied with args.parents) we'll prioritize
        // populations with conflict_ref
        include_populations = if !include_populations.is_empty() {
            debug!("Prioritizing parents: {:?}", &include_populations);
            include_populations
        }
        // prioritize populations that have min_subs conflict_alt
        // ie. help resolve the conflict_alt
        else {
            // count up the number of conflict_alt by population
            let mut conflict_alt_count = BTreeMap::new();
            for sub in conflict_alt {
                if dataset.mutations.contains_key(sub) {
                    let populations = dataset.mutations[sub]
                        .iter()
                        .filter(|pop| !exclude_populations.contains(pop))
                        .collect_vec();
                    for pop in populations {
                        *conflict_alt_count.entry(pop).or_insert(0) += 1
                    }
                }
            }

            let conflict_alt_populations = conflict_alt_count
                .into_iter()
                .filter(|(_pop, count)| *count >= args.min_subs)
                .map(|(pop, _count)| pop)
                .cloned()
                .collect_vec();

            // trunclate list for display
            let display_populations = if conflict_alt_populations.len() <= 10 {
                conflict_alt_populations.iter().join(", ")
            } else {
                format!("{} ...", conflict_alt_populations[0..10].iter().join(", "),)
            };

            debug!("Prioritizing conflict_alt resolution: {display_populations:?}");
            conflict_alt_populations
        };

        // --------------------------------------------------------------------
        // Search Dataset #1 (Full Coordinate Range)
        // --------------------------------------------------------------------

        // find the next parent candidate based on the full coordinate range
        let coord_min = coordinates.iter().min().unwrap();
        let coord_max = coordinates.iter().max().unwrap();
        let coord_range = (coord_min.to_owned()..coord_max.to_owned()).collect_vec();

        debug!("Searching based on coordinate range: {coord_min} - {coord_max}");

        let parent_candidate =
            dataset.search(sequence, Some(&include_populations), Some(&coord_range));

        // if the search found parents, check for recombination
        if let Ok(parent_candidate) = parent_candidate {
            // remove this parent from future searches
            exclude_populations.push(parent_candidate.consensus_population.clone());
            //designated_parents.retain(|p| p != &parent_candidate.consensus_population);

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
            if let Ok(detect_result) = detect_result {
                num_parents += 1;
                // reset the iter counter
                num_iter = 0;
                parents.push(parent_candidate);
                recombination = detect_result;
                continue;
            }
        }

        // --------------------------------------------------------------------
        // Search Dataset #2 (Precise Coordinates)

        debug!("Searching based on precise coordinates: {coordinates:?}");

        let parent_candidate =
            dataset.search(sequence, Some(&include_populations), Some(&coordinates));

        // if the search found parents, check for recombination
        if let Ok(parent_candidate) = parent_candidate {
            // remove this parent from future searches
            exclude_populations.push(parent_candidate.consensus_population.clone());
            //designated_parents.retain(|p| p != &parent_candidate.consensus_population);

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
            if let Ok(detect_result) = detect_result {
                num_parents += 1;
                // reset the iter counter
                num_iter = 0;
                parents.push(parent_candidate);
                recombination = detect_result;
                continue;
            }
        }
    }
}