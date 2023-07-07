use crate::cli::RunArgs;
use crate::dataset::{Dataset, SearchResult};
use crate::recombination::{detect_recombination, Recombination};
use crate::sequence::Sequence;
use color_eyre::eyre::{Report, Result};
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
    // initialize the two return values
    let mut parents = Vec::new();
    let mut recombination = Recombination::new(args.sequence);

    if args.max_parents == 0 {
        return Ok((parents, recombination));
    }

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
    (parents, recombination) = search_secondary(parents, &args)?;

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
    let mut exclude_populations = Vec::new();

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
            debug!(
                "Maximum number of parents reached ({}), stopping parent search.",
                args.max_parents
            );
            break;
        }
        if num_iter >= args.max_iter {
            debug!("Maximum number of iterations reached ({num_iter}), stopping parent search.");
            break;
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
        // EXCLUDE descendants/ancestors(?) of known parents

        for parent in &parents {
            // exclude descendants
            let descendants = args
                .dataset
                .phylogeny
                .get_descendants(&parent.consensus_population)?;
            let descendants_to_exclude = descendants
                .into_iter()
                .filter(|pop| !exclude_populations.contains(pop))
                .collect::<Vec<_>>();
            exclude_populations.extend(descendants_to_exclude);

            // exclude ancestors
            let ancestors = args
                .dataset
                .phylogeny
                .get_ancestors(&parent.consensus_population)?;
            // because of possible recombination, ancestors is a vector of vectors
            // to reflect multiple parents and paths to the root.
            // just flatten them for all our purposes here
            let ancestors_to_exclude = ancestors
                .into_iter()
                .flatten()
                .filter(|pop| !exclude_populations.contains(pop))
                .collect::<Vec<_>>();
            exclude_populations.extend(ancestors_to_exclude);

            // filter designated parents to prioritize search
            designated_parents = designated_parents
                .into_iter()
                .filter(|p| *p != parent.consensus_population)
                .collect_vec();
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
            .map(|(pop, _count)| pop.to_owned())
            .collect::<Vec<_>>();
        exclude_populations.extend(populations_to_exclude);

        // --------------------------------------------------------------------
        // Include populations that help resolve the conflict_alt

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

        // keep any populations that have min_subs conflict_alt
        let include_populations = conflict_alt_count
            .into_iter()
            .filter(|(_pop, count)| *count >= args.min_subs)
            .map(|(pop, _count)| pop)
            .cloned()
            .collect::<Vec<_>>();

        // --------------------------------------------------------------------
        // Search Dataset #1 (Full Coordinate Range)

        // find the next parent candidate based on the full coordinate range
        let coordinates_min = coordinates.iter().min().unwrap();
        let coordinates_max = coordinates.iter().max().unwrap();
        let coordinates_range =
            (coordinates_min.to_owned()..coordinates_max.to_owned()).collect_vec();

        // prioritize known/designated parents first
        let search_result = if !designated_parents.is_empty() {
            debug!("Prioritizing remaining designated parents: {designated_parents:?}");
            args.dataset.search(
                args.sequence,
                Some(&designated_parents),
                Some(&coordinates_range),
            )?
        }
        // otherwise prioritize populations that solve conflicts
        else {
            args.dataset.search(
                args.sequence,
                Some(&include_populations),
                Some(&coordinates_range),
            )?
        };
        recombination = detect_recombination(&parents, Some(&search_result), args)?;

        // if the recombination search succeeded,
        if !recombination.breakpoints.is_empty() {
            num_parents += 1;
            parents.push(search_result);
            continue;
        }

        // --------------------------------------------------------------------
        // Search Dataset #2 (Specific Listed Coordinates)

        // prioritize known/designated parents first
        let search_result = if !designated_parents.is_empty() {
            debug!("Prioritizing remaining designated parents: {designated_parents:?}");
            args.dataset.search(
                args.sequence,
                Some(&designated_parents),
                Some(&coordinates),
            )?
        }
        // otherwise prioritize populations that solve conflicts
        else {
            args.dataset.search(
                args.sequence,
                Some(&include_populations),
                Some(&coordinates),
            )?
        };
        recombination = detect_recombination(&parents, Some(&search_result), args)?;

        // if the recombination search succeeded,
        if !recombination.breakpoints.is_empty() {
            num_parents += 1;
            parents.push(search_result);
            continue;
        }
        // if failed, exclude search_result top_populations from next iteration?
        // do we also want to exclude failed search #1 top_populations?
        else {
            exclude_populations.extend(search_result.top_populations);
        }
    }

    Ok((parents, recombination))
}
