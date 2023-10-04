use crate::cli::run;
use crate::dataset::{Dataset, SearchResult};
use crate::recombination::{detect_recombination, Recombination};
use crate::sequence::Sequence;
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;
use log::{debug, warn};
use std::collections::BTreeMap;

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

// Search for all parents (primary and secondary parents).
pub fn all_parents<'seq>(
    sequence: &'seq Sequence,
    dataset: &Dataset,
    best_match: &mut SearchResult,
    args: &run::Args,
) -> Result<(Vec<SearchResult>, Recombination<'seq>), Report> {

    // This is not a good place to check this argument, should be further upstream?
    if args.max_parents == 0 {
        return Err(eyre!("Parameter max_parents is set to 0."));
    }

    let mut edge_case = false;
    let mut fallback_bias = true;

    // ----------------------------------------------------------------------------
    // Naive Search
    // ----------------------------------------------------------------------------

    // use best match as the first parent
    let parents = vec![best_match.clone()];
    // try to search for secondary parents
    let mut result = secondary_parents(sequence, dataset, best_match, &parents, &args);
    // examine the result, we won't run the fallback bias search if multiple parents were found
    if let Ok((parents, _)) = &result {
        if parents.len() > 1 {
            fallback_bias = false
        }
    }

    // ----------------------------------------------------------------------------
    // Biased Search
    // ----------------------------------------------------------------------------

    if !args.naive && fallback_bias {

        warn!("Naive recombination search failed, trying biased search.");

        let mut exclude_populations: Vec<String> = vec![];
        let mut include_populations: Vec<String> = vec![];
        let mut parents = vec![];
        let mut args = args.clone();

        // Bias #1. Use phylogenetic information about designated parents
        if let Some(recombinant) = &best_match.recombinant {
            let designated_parents = dataset.phylogeny.get_parents(recombinant)?;
            debug!("Applying bias #1, prioritizing designated parents: {designated_parents:?}");
            include_populations = designated_parents;
        }

        // Bias #2. Use edge case parameters (if applicable)
        let edge_case_args = dataset
            .edge_cases
            .iter()
            .find(|e| e.population == Some(best_match.consensus_population.to_string()));

        if let Some(edge_case_args) = edge_case_args {
            debug!("Applying Search Bias #2, edge case found: {edge_case_args:?}");
            edge_case = true;
            args = args.apply_edge_case(edge_case_args)?;
            // If parents are part of the edge case parameter, supply them
            if let Some(parents) = &edge_case_args.parents {
                include_populations = parents.clone();
            }
        }        

        // Bias #3. Use parents from CLI Args
        if let Some(parents_cli) = &args.parents {
            debug!("Applying Search Bias #3, parents from CLI: {parents_cli:?}");
            include_populations = parents_cli.clone();
        }

        // ----------------------------------------------------------------------------
        // Primary parent

        // Possibly different parameters for populations
        let parent_primary = if include_populations.is_empty() {
            dataset.search(sequence, None, None)?
        } else {
            dataset.search(sequence, Some(&include_populations), None)?
        };
        parents.push(parent_primary);

        // ----------------------------------------------------------------------------
        // Secondary parents ( 2 : max_parents)
        result = secondary_parents(sequence, dataset, best_match, &parents, &args);
    }

    let (parents, mut recombination) = result?;
    recombination.edge_case = edge_case;

    // check if observed parents match expected parents
    // ex. XBL (before designated) was XBB ( BA.2.75 and XBB.1.5.57)
    // which conflicts with expected XBB parents BJ.1 (BA.2.10.1) and BM.1.1.1 (BA.2.75.3)
    // todo!() decide on the order, should descendants be from observed or expected?
    if let Some(recombinant) = &best_match.recombinant {

        let expected_parents = dataset.phylogeny.get_parents(&recombinant)?;
        let mut observed_parents_descendants = recombination.parents
            .iter()
            .flat_map(|p| dataset.phylogeny.get_descendants(p).unwrap())
            .unique()
            .collect_vec();

        observed_parents_descendants.retain(|p| expected_parents.contains(p));
        if observed_parents_descendants.len() != expected_parents.len(){
            best_match.recombinant = Some("novel".to_string());
        }
    }

    Ok((parents, recombination))
}

// Search for the secondary recombination parent(s).
pub fn secondary_parents<'seq>(
    sequence: &'seq Sequence,
    dataset: &Dataset,
    best_match: &SearchResult,
    parents: &Vec<SearchResult>,
    args: &run::Args,
) -> Result<(Vec<SearchResult>, Recombination<'seq>), Report> {
    let mut recombination = Recombination::new(sequence);
    let mut parents = parents.clone();
    let mut num_iter = 0;
    let mut num_parents = 1;

    let mut exclude_populations: Vec<String> = Vec::new();
    let mut designated_parents = Vec::new();

    // // running a biased search, prioritize designated parents
    // if !args.naive {
    //     if let Some(recombinant) = &best_match.recombinant {
    //         designated_parents = dataset.phylogeny.get_parents(recombinant)?;
    //         // exclude self and descendants from secondary parent
    //         // ex. XR might try to test against itself as a parent
    //         let descendants = dataset.phylogeny.get_descendants(recombinant)?;
    //         exclude_populations.extend(descendants)
    //     }
    // }

    loop {
        // --------------------------------------------------------------------
        // Loop Break Check: Simple
        //
        // Check if we can break out of the loop based on simple checks like the
        // max number of iterations max number of parents achieved.

        if num_parents >= args.max_parents {
            debug!("Maximum parents reached ({num_parents}).");
            return Ok((parents, recombination));
        }
        if num_iter >= args.max_iter {
            debug!("Maximum iterations reached ({num_iter}).");

            // If we found parents, fine to max out iter
            if num_parents > 1 {
                return Ok((parents, recombination));
            }
            // If we didn't find parents, this is a failure
            else {
                return Err(eyre!("No secondary parents found."));
            }
        }
        num_iter += 1;

        debug!("parent_{}: iteration {num_iter}", num_parents + 1);

        // --------------------------------------------------------------------
        // Conflict Checks
        //
        // Check if we have resolved all mutation conflicts between the query
        // sequence and recombination parents found so far.

        // identify all substitutions found in all parents so far
        let parent_substitutions = parents
            .iter()
            .flat_map(|parent| &parent.substitutions)
            .unique()
            .collect_vec();

        // identify conflict_alt (alt base) in sequence not resolved by any parent
        let conflict_alt = sequence
            .substitutions
            .iter()
            .filter(|sub| !parent_substitutions.contains(sub))
            .collect_vec();

        debug!("conflict_alt: {}", conflict_alt.iter().join(", "));

        // identify conflict_ref (ref base) in sequence not resolved by any parent
        // first collect all the conflict_ref in all parents
        let mut conflict_ref = parents
            .iter()
            .flat_map(|p| &p.conflict_ref[&p.consensus_population])
            .unique()
            .collect_vec();

        // next, check for any conflict_ref that is NOT resolved by another parent
        // a conflict_ref is resolved if one parent DOES NOT have it in subs
        // this is not ideal, because it does not account for missing and deletions
        let mut conflict_ref_unresolved = Vec::new();
        for sub in conflict_ref {
            let mut is_resolved = false;
            for parent in &parents {
                if !parent.substitutions.contains(sub) {
                    is_resolved = true;
                }
            }
            if !is_resolved {
                conflict_ref_unresolved.push(sub);
            }
        }

        // filter the conflict_ref down to just unresolved
        conflict_ref = conflict_ref_unresolved;
        debug!("conflict_ref: {}", conflict_ref.iter().join(", "));

        // --------------------------------------------------------------------
        // Loop Break Check: Conflict

        // if conflict_ref.is_empty() {
        //     debug!("Sufficient conflict_ref resolution reached, stopping parent search.");
        //     break;
        // }
        if conflict_alt.len() < args.min_subs {
            debug!("Sufficient conflict_alt resolution reached, stopping parent search.");
            // If we found parents, fine to max out iter
            if num_parents > 1 {
                return Ok((parents, recombination));
            } else {
                return Err(eyre!("No secondary parents found."));             
            }
        }

        // Collect the conflict/unresolved coordinates
        let conflict_ref_coord = conflict_ref
            .iter()
            .map(|sub| &sub.coord)
            .cloned()
            .collect_vec();
        let conflict_alt_coord = conflict_alt
            .iter()
            .map(|sub| &sub.coord)
            .cloned()
            .collect_vec();
        let mut coordinates = conflict_ref_coord;
        coordinates.extend(conflict_alt_coord);
        coordinates.sort();
        debug!("coordinates: {coordinates:?}");

        // filter designated parents to prioritize search
        for parent in &parents {
            designated_parents.retain(|p| p != &parent.consensus_population);
        }

        // --------------------------------------------------------------------
        // Exclude populations that have substitutions at ALL of the conflict_ref

        let mut conflict_ref_count = BTreeMap::new();
        for sub in conflict_ref.iter() {
            if dataset.mutations.contains_key(sub) {
                let populations = dataset.mutations[sub]
                    .iter()
                    .filter(|pop| !exclude_populations.contains(pop))
                    .collect_vec();
                for pop in populations {
                    *conflict_ref_count.entry(pop).or_insert(0) += 1
                }
            }
        }
        let populations_to_exclude = conflict_ref_count
            .into_iter()
            .filter(|(_pop, count)| *count == conflict_ref.len())
            .map(|(pop, _count)| pop.clone())
            .collect_vec();
        exclude_populations.extend(populations_to_exclude);

        // --------------------------------------------------------------------
        // Include Populations

        // prioritize known/designated parents first in search
        let include_populations = if !designated_parents.is_empty() {
            debug!("Prioritizing designated parents: {designated_parents:?}");
            designated_parents.clone()
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

            //debug!("Prioritizing conflict_alt resolution: {conflict_alt_count:?}");
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
            designated_parents.retain(|p| p != &parent_candidate.consensus_population);

            // check for recombination
            let detect_result = detect_recombination(
                sequence,
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
            designated_parents.retain(|p| p != &parent_candidate.consensus_population);

            // check for recombination
            let detect_result = detect_recombination(
                sequence,
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
            }
        }
    }
}
