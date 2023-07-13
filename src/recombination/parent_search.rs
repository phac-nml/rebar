use crate::cli::RunArgs;
use crate::dataset::{Dataset, SearchResult};
use crate::recombination::{detect_recombination, Recombination};
use crate::sequence::Sequence;
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;
use log::debug;
use std::collections::BTreeMap;

// ----------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------

pub struct Args<'data, 'seq, 'best> {
    pub dataset: &'data Dataset,
    pub sequence: &'seq Sequence,
    pub best_match: &'best SearchResult,
    pub max_parents: usize,
    pub max_iter: usize,
    pub min_consecutive: usize,
    pub min_length: usize,
    pub min_subs: usize,
}

impl<'data, 'seq, 'best, 'args> Args<'data, 'seq, 'best> {
    pub fn new(
        dataset: &'data Dataset,
        sequence: &'seq Sequence,
        best_match: &'best SearchResult,
        args: &'args RunArgs,
    ) -> Self {
        Args {
            dataset,
            sequence,
            best_match,
            max_parents: args.max_parents,
            max_iter: args.max_iter,
            min_consecutive: args.min_consecutive,
            min_length: args.min_length,
            min_subs: args.min_subs,
        }
    }
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

// Search for all parents (primary and secondary recombination).
pub fn search_all<'seq>(
    args: Args<'_, 'seq, '_>,
) -> Result<(Vec<SearchResult>, Recombination<'seq>), Report> {
    if args.max_parents == 0 {
        return Err(eyre!("Parameter max_parents is set to 0."));
    }

    // initialize parents to return
    let mut parents = Vec::new();

    // don't restrict which populations are excluded/included
    // by default, dataset.search will look through all populations
    let mut exclude_populations = Vec::new();
    let mut include_populations = Vec::new();

    // For known recombinants, prioritize designated parent
    let recombinant = &args.best_match.recombinant;
    if let Some(recombinant) = recombinant {
        include_populations = args.dataset.phylogeny.get_parents(recombinant)?;
    }

    // Primary parent
    let parent_primary =
        search_primary(&args, &mut include_populations, &mut exclude_populations)?;
    parents.push(parent_primary);

    // Secondary parents ( 2 : max_parents)
    // this function consumes `parents`, modifies it, then returns it
    let (parents, recombination) = search_secondary(parents, &args)?;

    Ok((parents, recombination))
}

// Search for the primary parent.
pub fn search_primary(
    args: &Args,
    include_populations: &mut Vec<String>,
    exclude_populations: &mut Vec<String>,
) -> Result<SearchResult, Report> {
    debug!("Searching for Parent 1");

    // If this is a known recombinant, exclude self and the recombinant's descendants from parent search.
    let recombinant = args.best_match.recombinant.clone();
    let search_result = if let Some(recombinant) = recombinant {
        // exclude self
        exclude_populations.push(args.best_match.consensus_population.clone());
        // exclude descendants
        let descendants = args.dataset.phylogeny.get_descendants(&recombinant)?;
        exclude_populations.extend(descendants);
        *include_populations = include_populations
            .iter()
            .filter(|pop| !exclude_populations.contains(pop))
            .map(|pop| pop.to_owned())
            .collect::<Vec<_>>();
        args.dataset
            .search(args.sequence, Some(include_populations), None)?
    }
    // Not a known recombinant, just use best_match/consensus as parent 1
    else {
        args.best_match.to_owned()
    };

    Ok(search_result)
}

// Search for the secondary recombination parent(s).
pub fn search_secondary<'seq>(
    mut parents: Vec<SearchResult>,
    args: &Args<'_, 'seq, '_>,
) -> Result<(Vec<SearchResult>, Recombination<'seq>), Report> {
    let mut recombination = Recombination::new(args.sequence);
    let mut num_iter = 0;
    let mut num_parents = 1;

    let mut exclude_populations: Vec<String> = Vec::new();

    // For known recombinants, prioritize designated parents
    let mut designated_parents = Vec::new();
    let recombinant = &args.best_match.recombinant;
    if let Some(recombinant) = recombinant {
        designated_parents = args.dataset.phylogeny.get_parents(recombinant)?;
    }

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
            // If we found parents, fine to max out iter
            if num_parents > 1 {
                debug!("Maximum iterations reached ({num_iter}).");
                return Ok((parents, recombination));
            }
            // If we didn't find parents, this is a failure
            else {
                return Err(eyre!("Maximum iterations reached ({num_iter}) with no secondary parent found."));
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
        let conflict_alt = args
            .sequence
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
            break;
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

        // --------------------------------------------------------------------
        // Exclude Populations

        for parent in &parents {
            // // exclude descendants of known parents
            // let descendants = args
            //     .dataset
            //     .phylogeny
            //     .get_descendants(&parent.consensus_population)?;
            // let descendants_to_exclude = descendants
            //     .into_iter()
            //     .filter(|pop| !exclude_populations.contains(pop))
            //     .collect::<Vec<_>>();
            // exclude_populations.extend(descendants_to_exclude);

            // // exclude ancestors of known parents
            // let ancestors = args
            //     .dataset
            //     .phylogeny
            //     .get_ancestors(&parent.consensus_population)?;
            // // because of possible recombination, ancestors is a vector of vectors
            // // to reflect multiple parents and paths to the root.
            // // just flatten them for all our purposes here
            // let ancestors_to_exclude = ancestors
            //     .into_iter()
            //     .flatten()
            //     .filter(|pop| !exclude_populations.contains(pop))
            //     .collect::<Vec<_>>();
            // exclude_populations.extend(ancestors_to_exclude);

            // filter designated parents to prioritize search
            designated_parents.retain(|p| p != &parent.consensus_population);
        }

        // --------------------------------------------------------------------
        // Exclude populations that have substitutions at ALL of the conflict_ref

        let mut conflict_ref_count = BTreeMap::new();
        for sub in conflict_ref.iter() {
            if args.dataset.mutations.contains_key(sub) {
                let populations = args.dataset.mutations[sub]
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
                if args.dataset.mutations.contains_key(sub) {
                    let populations = args.dataset.mutations[sub]
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

            debug!("Prioritizing conflict_alt resolution: {conflict_alt_populations:?}");
            conflict_alt_populations
        };

        // --------------------------------------------------------------------
        // Search Dataset #1 (Full Coordinate Range)

        // find the next parent candidate based on the full coordinate range
        let coord_min = coordinates.iter().min().unwrap();
        let coord_max = coordinates.iter().max().unwrap();
        let coord_range = (coord_min.to_owned()..coord_max.to_owned()).collect_vec();

        debug!("Searching based on coordinate range: {coord_min} - {coord_max}");

        let search_result = args.dataset.search(
            args.sequence,
            Some(&include_populations),
            Some(&coord_range),
        );

        // if the search found parents, check for recombination
        if let Ok(search_result) = search_result {
            // remove this parent from future searches
            exclude_populations.push(search_result.consensus_population.clone());
            designated_parents.retain(|p| p != &search_result.consensus_population);

            // check for recombination
            let detect_result =
                detect_recombination(&parents, Some(&search_result), args);

            // if successful, add this parent to the list and update recombination
            if let Ok(detect_result) = detect_result {
                num_parents += 1;
                parents.push(search_result);
                recombination = detect_result;
                continue;
            }
        }

        // --------------------------------------------------------------------
        // Search Dataset #2 (Precise Coordinates)

        debug!("Searching based on precise coordinates: {coordinates:?}");

        let search_result = args.dataset.search(
            args.sequence,
            Some(&include_populations),
            Some(&coordinates),
        );

        // if the search found parents, check for recombination
        if let Ok(search_result) = search_result {
            // remove this parent from future searches
            exclude_populations.push(search_result.consensus_population.clone());
            designated_parents.retain(|p| p != &search_result.consensus_population);

            // check for recombination
            let detect_result =
                detect_recombination(&parents, Some(&search_result), args);
            // if successful, add this parent to the list and update recombination
            if let Ok(detect_result) = detect_result {
                num_parents += 1;
                parents.push(search_result);
                recombination = detect_result;
                continue;
            }
        }
    }

    Ok((parents, recombination))
}
