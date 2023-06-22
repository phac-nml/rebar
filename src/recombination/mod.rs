use crate::dataset::{Dataset, SearchResult};
use crate::cli::RunArgs;                         // object for run function params
use crate::sequence::{Sequence, Substitution};
use color_eyre::eyre::{Report, Result};
use itertools::Itertools;
use log::debug;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::default::Default;
use tabled::builder::Builder;
use tabled::settings::Style;

// ----------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Breakpoint

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Breakpoint {
    pub start: usize,
    pub end: usize,
}

// ----------------------------------------------------------------------------
// Direction

pub enum Direction {
    Forward,
    Reverse,
}

// ----------------------------------------------------------------------------
// Region

#[derive(Clone, Debug, Deserialize, Serialize)]
#[allow(dead_code)]
pub struct Region {
    pub start: usize,
    pub end: usize,
    pub origin: String,
    #[serde(skip_serializing)]
    pub substitutions: Vec<Substitution>,
}

// ----------------------------------------------------------------------------
// Recombination

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Recombination {
    pub breakpoints: Vec<Breakpoint>,
    pub regions: BTreeMap<usize, Region>,
    pub table: Vec<Vec<String>>,
}

impl Recombination {
    pub fn new() -> Self {
        Recombination {
            breakpoints: Vec::new(),
            regions: BTreeMap::new(),
            table: Vec::new(),
        }
    }
}

impl Default for Recombination {
    fn default() -> Self {
        Self::new()
    }
}

// ----------------------------------------------------------------------------
// Find Parents

pub struct ParentSearchArgs<'d, 's, 'b> {
    pub dataset: &'d Dataset,
    pub sequence: &'s Sequence,
    pub best_match: &'b SearchResult,
    pub max_parents: usize,
    pub max_iter: usize,
    pub min_consecutive: usize,
    pub min_length: usize,
    pub min_subs: usize,
}

impl<'d, 's, 'b, 'r> ParentSearchArgs<'d, 's, 'b> {
    pub fn new(dataset: &'d Dataset, sequence: &'s Sequence, best_match: &'b SearchResult, args: &'r RunArgs) -> Self {
        ParentSearchArgs {
            dataset,
            sequence,
            best_match,
            max_parents : args.max_parents,
            max_iter : args.max_iter,
            min_consecutive : args.min_consecutive,
            min_length : args.min_length,
            min_subs : args.min_subs,
        }
    }
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

// Search function specific to the secondary parent(s).
pub fn parent_search_secondary(
    sequence: &Sequence, 
    mut parents: Vec<SearchResult>, 
    args: &ParentSearchArgs
) -> Result<(Vec<SearchResult>, Recombination), Report> {

    let mut recombination = Recombination::new();

    let mut num_iter = 0;
    let mut num_parents = 1;
    let mut exclude_populations = Vec::new();

    loop {

        // --------------------------------------------------------------------
        // Loop Break Check: Simple
        //
        // Check if we can break out of the loop based on simple checks like the
        // max number of iterations max number of parents achieved.

        if num_parents >= args.max_parents {
            debug!("Maximum number of parents reached ({}), stopping parent search.", args.max_parents);
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
            .map(|parent| &parent.substitutions)
            .flatten()
            .unique()
            .collect_vec();

        // identify conflict_alt (alt base) in sequence not resolved by any parent
        let conflict_alt = sequence.substitutions
            .iter()
            .filter(|sub| !parent_substitutions.contains(&sub))
            .collect_vec();

        debug!("conflict_alt: {}", conflict_alt.iter().join(", "));

        // identify conflict_ref (ref base) in sequence not resolved by any parent
        // first collect all the conflict_ref in all parents
        let mut conflict_ref = parents
            .iter()
            .map(|p| &p.conflict_ref[&p.consensus_population])
            .flatten()
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
        let conflict_ref_coord = conflict_ref.iter().map(|sub| &sub.coord).cloned().collect_vec();
        let conflict_alt_coord = conflict_alt.iter().map(|sub| &sub.coord).cloned().collect_vec();
        let mut coordinates = conflict_ref_coord;
        coordinates.extend(conflict_alt_coord);
        coordinates.sort();
        debug!("coordinates: {coordinates:?}");

        // --------------------------------------------------------------------
        // EXCLUDE descendants/ancestors(?) of known parents

        for parent in &parents {

            // exclude descendants
            let descendants = args.dataset
                .phylogeny
                .get_descendants(&parent.consensus_population)?;
            let descendants_to_exclude = descendants
                .into_iter()
                .filter(|pop| !exclude_populations.contains(pop))
                .collect::<Vec<_>>();
            exclude_populations.extend(descendants_to_exclude); 

            // exclude ancestors
            let ancestors = args.dataset
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
        let coordinates_range = (coordinates_min.to_owned()..coordinates_max.to_owned()).collect_vec();

        let search_result = args.dataset.search(sequence, Some(&include_populations), Some(&coordinates_range))?;        
        recombination = detect_recombination(&sequence, &parents, Some(&search_result), &args)?;

        // if the recombination search succeeded, 
        if !recombination.breakpoints.is_empty() {
            num_parents += 1;
            parents.push(search_result);
            continue;
        } 

        // --------------------------------------------------------------------
        // Search Dataset #2 (Specific Listed Coordinates)

        let search_result = args.dataset.search(sequence, Some(&include_populations), Some(&coordinates))?;        
        recombination = detect_recombination(&sequence, &parents, Some(&search_result), &args)?;        
    
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

pub fn parent_search(args: ParentSearchArgs) -> Result<(Vec<SearchResult>, Recombination), Report> {

    // initialize the two return values
    let mut parents = Vec::new();
    let mut recombination = Recombination::new();

    if args.max_parents == 0 {
        return Ok((parents, recombination));
    }

    // by default, include all dataset populations
    let mut include_populations = args.dataset
        .populations
        .keys()
        .map(|pop| pop.to_owned())
        .collect::<Vec<_>>();

    // by default, don't exclude any populations
    let mut exclude_populations = Vec::new();

    // Primary parent
    let parent_primary = parent_search_primary(&args, &mut include_populations, &mut exclude_populations)?;
    parents.push(parent_primary);

    // Secondary parents ( 2 : max_parents)
    // this function consumes `parents`, modifies it, then returns it
    (parents, recombination) = parent_search_secondary(&args.sequence, parents, &args)?;

    Ok((parents, recombination))
}

pub fn detect_recombination(
    sequence: &Sequence,
    parents: &Vec<SearchResult>,
    search_result: Option<&SearchResult>,
    args: &ParentSearchArgs,
) -> Result<Recombination, Report> {
    let mut recombination = Recombination::new();

    // if no search_result was provided, just use first parent
    let search_result_default = &parents[0];
    let search_result = match search_result {
        Some(summary) => summary,
        None => search_result_default,
    };

    // Create a table where rows are coordinates and columns are
    // coord, parent, Reference, <parents...>, <search_result> <sequence>

    let mut table_rows = Vec::new();

    // Identify which subs are non-bi-allelic, these will wind up being
    // duplicate rows, which we'll need to reconcile and collapse
    let mut all_subs = Vec::new();
    for parent in parents {
        all_subs.extend(parent.substitutions.to_owned());
    }
    all_subs.extend(search_result.substitutions.to_owned());
    all_subs.sort();

    // --------------------------------------------------------------------
    // Identify Substitution Parental Origins
    // --------------------------------------------------------------------

    let all_coords = all_subs.iter().map(|sub| sub.coord).collect::<Vec<_>>();
    // Deduplicate coords (might be non-bi-allelic)
    let coords = all_coords.into_iter().unique().collect::<Vec<_>>();
    let mut privates = Vec::new();

    for coord in coords {
        // add coord as first column
        let mut row = vec![coord.to_string()];

        // init base origins (could be multiple)
        let mut origins = Vec::new();

        // add reference base as second column
        let ref_base = all_subs
            .iter()
            .filter(|sub| sub.coord == coord)
            .map(|sub| sub.reference)
            .next()
            .unwrap();
        row.push(ref_base.to_string());

        // Base in Sample
        let sample_base = sequence
            .substitutions
            .iter()
            .filter(|sub| sub.coord == coord)
            .map(|sub| sub.alt)
            .next()
            .unwrap_or(ref_base);

        // Add the known parents as next columns
        let mut parent_bases = Vec::new();

        for parent in parents {
            let parent_base = parent
                .substitutions
                .iter()
                .filter(|sub| sub.coord == coord)
                .map(|sub| sub.alt)
                .next()
                .unwrap_or(ref_base);

            if parent_base == sample_base {
                origins.push(parent.consensus_population.clone());
            }

            row.push(parent_base.to_string());
            parent_bases.push(parent_base);
        }

        // Add the current parent match that is being evaluated
        let match_base = search_result
            .substitutions
            .iter()
            .filter(|sub| sub.coord == coord)
            .map(|sub| sub.alt)
            .next()
            .unwrap_or(ref_base);

        if match_base == sample_base {
            origins.push(search_result.consensus_population.clone());
        }
        row.push(match_base.to_string());
        parent_bases.push(match_base);

        // Add the sample base
        row.push(sample_base.to_string());

        // Is this coord a discriminating site?
        // Remove subs that are identical between all parents
        parent_bases = parent_bases.into_iter().unique().collect();

        // If only 1 unique parent base was found, non-discriminating
        if parent_bases.len() == 1 {
            continue;
        }

        // If no origins were found, this is private
        if origins.is_empty() {
            let private = sequence.substitutions.iter().find(|sub| sub.coord == coord);
            privates.push(private);
        } else if origins.len() == 1 {
            row.push(origins[0].clone());
            table_rows.push(row);
        }
    }

    // --------------------------------------------------------------------
    // Group Substitutions into Parental Regions
    // --------------------------------------------------------------------

    // First: 5' -> 3', filter separately on min_consecutive then min_length
    let mut regions_5p = identify_regions(&table_rows)?;
    regions_5p =
        filter_regions(&regions_5p, Direction::Forward, args.min_consecutive, 0_usize)?;
    regions_5p = filter_regions(&regions_5p, Direction::Forward, 0_usize, args.min_length)?;
    debug!(
        "regions_5p: {}",
        serde_json::to_string(&regions_5p).unwrap()
    );

    // Second: 3' -> 5', filter separately on min_consecutive then min_length
    let mut regions_3p = identify_regions(&table_rows)?;
    regions_3p =
        filter_regions(&regions_3p, Direction::Reverse, args.min_consecutive, 0_usize)?;
    regions_3p = filter_regions(&regions_3p, Direction::Reverse, 0_usize, args.min_length)?;
    debug!(
        "regions_3p: {}",
        serde_json::to_string(&regions_3p).unwrap()
    );

    // Take the intersect of the 5' and 3' regions (ie. where they both agree)
    let mut regions_intersect = intersect_regions(&regions_5p, &regions_3p)?;
    // During intersection, it's possible that a region from 1 single parent
    // got broken up into multiple adjacent sections. Put it through the
    // filter again to collapse it. Direction doesn't matter now
    regions_intersect = filter_regions(
        &regions_intersect,
        Direction::Forward,
        args.min_consecutive,
        args.min_length,
    )?;
    debug!(
        "regions_intersect: {}",
        serde_json::to_string(&regions_intersect).unwrap()
    );

    // Make sure all the prev_parents + search_result have at least 1 region
    let region_origins = regions_intersect
        .values()
        .map(|region| region.origin.to_owned())
        .collect_vec();

    for parent in parents {
        if !region_origins.contains(&parent.consensus_population) {
            debug!(
                "No recombination detected for parent {}.",
                &parent.consensus_population
            );
            return Ok(recombination);
        }
    }
    if !region_origins.contains(&search_result.consensus_population) {
        debug!(
            "No recombination detected for parent {}.",
            &search_result.consensus_population
        );
        return Ok(recombination);
    }

    // --------------------------------------------------------------------
    // Minimum Substitutions Filter
    // --------------------------------------------------------------------

    // Check for the minimum number of subs from each parent
    let mut uniq_subs_count: BTreeMap<String, usize> = BTreeMap::new();
    // first count up uniq subs by parent
    for region in regions_intersect.values() {
        let origin = &region.origin;
        uniq_subs_count.entry(origin.clone()).or_insert(0);
        for sub in &region.substitutions {
            if sub.reference == sub.alt {
                continue;
            }
            *uniq_subs_count.entry(origin.clone()).or_default() += 1;
        }
    }
    debug!("uniq_subs_count: {uniq_subs_count:?}");
    // check if any parent failed the filter
    let mut min_subs_fail = false;
    for (parent, count) in uniq_subs_count {
        if count < args.min_subs {
            debug!(
                "Parent {} subs ({count}) do not meet the min_subs filter ({}).",
                parent,
                args.min_subs,
            );
            min_subs_fail = true;
        }
    }
    if min_subs_fail {
        debug!("No recombination detected, min_subs filter was not satisfied by all parents.");
        return Ok(recombination);
    }

    // --------------------------------------------------------------------
    // Debugging Table (post-filter)
    // --------------------------------------------------------------------

    // table headers
    let mut headers = vec!["coord".to_string(), "Reference".to_string()];
    for parent in parents {
        headers.push(parent.consensus_population.to_owned());
    }
    headers.push(search_result.consensus_population.to_owned());
    headers.push(sequence.id.to_owned());
    headers.push("origin".to_string());

    // filter the table for the region subs
    let mut region_sub_coords = Vec::new();
    for region in regions_intersect.values() {
        let sub_coords = region
            .substitutions
            .iter()
            .map(|sub| sub.coord)
            .collect::<Vec<_>>();
        region_sub_coords.extend(sub_coords);
    }

    let mut table_rows_filter = Vec::new();
    for row in table_rows {
        let coord = row[0].parse::<usize>().unwrap();
        if region_sub_coords.contains(&coord) {
            table_rows_filter.push(row);
        }
    }

    let mut table_builder = Builder::default();
    table_builder.set_header(headers.clone());

    // table rows
    for row in &table_rows_filter {
        table_builder.push_record(row.to_owned());
    }

    let mut table = table_builder.build();

    // pretty print table for debugging
    debug!(
        "Recombination table:\n{}",
        table.with(Style::sharp()).to_string()
    );

    // --------------------------------------------------------------------
    // Breakpoints
    // --------------------------------------------------------------------

    let breakpoints = identify_breakpoints(&regions_intersect)?;
    debug!(
        "breakpoints: {}",
        serde_json::to_string(&breakpoints).unwrap()
    );

    // --------------------------------------------------------------------
    // Update
    // --------------------------------------------------------------------

    // for the table_rows, combine header and data
    let mut table = vec![headers];
    table.extend(table_rows_filter);

    // update all the attributes
    recombination.regions = regions_intersect;
    recombination.breakpoints = breakpoints;
    recombination.table = table;

    Ok(recombination)
}

// Search function specific to the primary parent.
pub fn parent_search_primary(
    args: &ParentSearchArgs, 
    include_populations: &mut Vec<String>,
    exclude_populations: &mut Vec<String>,
) -> Result<SearchResult, Report> {

    debug!("Searching for Parent 1");

    // If this is a known recombinant, exclude self the recombinant's descendants from parent search.
    let recombinant = args.best_match.recombinant.clone();
    let search_result = if let Some(recombinant) = recombinant {
        // exclude self
        exclude_populations.push(args.best_match.consensus_population.clone());
        // exclude descendants
        let descendants = args.dataset.phylogeny.get_descendants(&recombinant)?;
        exclude_populations.extend(descendants);
        *include_populations = include_populations
            .iter()
            .filter(|pop| !exclude_populations.contains(&pop))
            .map(|pop| pop.to_owned())
            .collect::<Vec<_>>();
        args.dataset.search(args.sequence, Some(include_populations), None)?
    }
    // Not a known recombinant, just use best_match/consensus as parent 1
    else {
        args.best_match.to_owned()
    };

    Ok(search_result)
}

pub fn identify_regions(
    table_rows: &Vec<Vec<String>>,
) -> Result<BTreeMap<usize, Region>, Report> {
    let mut origin_prev: Option<String> = None;
    let mut regions = BTreeMap::new();
    let mut start = 0;

    for row in table_rows {
        let coord = row[0].parse::<usize>().unwrap();
        let origin = row[row.len() - 1].to_string();
        let reference = row[1].chars().next().unwrap();
        let alt = row[row.len() - 2].chars().next().unwrap();
        let substitutions = vec![Substitution {
            coord,
            reference,
            alt,
        }];

        // start of new region, either first or origin changes
        if origin_prev.is_none() || origin_prev != Some(origin.clone()) {
            start = coord;
            let region = Region {
                start,
                end: coord,
                origin: origin.clone(),
                substitutions,
            };
            regions.insert(start, region);
        }
        // same origin, region continues. update end and subs
        else {
            let region = regions.get_mut(&start).unwrap();
            region.end = coord;
            region.substitutions.extend(substitutions)
        }

        origin_prev = Some(origin);
    }

    Ok(regions)
}

/// Filter recombinant regions based on the length and consecutive bases.
pub fn filter_regions(
    regions: &BTreeMap<usize, Region>,
    direction: Direction,
    min_consecutive: usize,
    min_length: usize,
) -> Result<BTreeMap<usize, Region>, Report> {
    let mut regions_filter = BTreeMap::new();
    let mut origin_prev: Option<String> = None;
    let mut start_prev: Option<usize> = None;

    let start_coords = match direction {
        Direction::Forward => regions.keys().collect::<Vec<_>>(),
        Direction::Reverse => regions.keys().rev().collect::<Vec<_>>(),
    };

    for start in start_coords {
        let region = regions.get(start).unwrap();
        let num_consecutive = region.substitutions.len();
        let region_length = (region.end - region.start) + 1;

        // start of new region, either first or origin changes
        if origin_prev.is_none() || origin_prev != Some(region.origin.clone()) {
            // is the new parental region long enough?
            if num_consecutive >= min_consecutive && region_length >= min_length {
                regions_filter.insert(region.start, region.to_owned());
                origin_prev = Some(region.origin.clone());
                start_prev = Some(region.start);
            }
        }
        // same origin, region continues. update subs and start or end
        else {
            match direction {
                // when going forward, we update the region
                Direction::Forward => {
                    if let Some(start_prev) = start_prev {
                        let region_update = regions_filter.get_mut(&start_prev).unwrap();
                        region_update
                            .substitutions
                            .extend(region.substitutions.clone());
                        region_update.end = region.end;
                    }
                }
                // when going backward, we remove and replace regions
                Direction::Reverse => {
                    if let Some(start_prev) = start_prev {
                        let mut region_new =
                            regions_filter.get(&start_prev).unwrap().to_owned();
                        region_new
                            .substitutions
                            .extend(region.substitutions.clone());
                        region_new.substitutions.sort();
                        region_new.start = region.start;

                        // remove old region from filtered map
                        regions_filter.remove(&start_prev);
                        // add new region with updated start coord
                        regions_filter.insert(region.start, region_new);
                    }

                    // for reverse, update new start position
                    start_prev = Some(region.start);
                }
            }
        }
    }

    Ok(regions_filter)
}

/// Find the intersect between two regions.
pub fn intersect_regions(
    regions_1: &BTreeMap<usize, Region>,
    regions_2: &BTreeMap<usize, Region>,
) -> Result<BTreeMap<usize, Region>, Report> {
    let mut regions_intersect = BTreeMap::new();

    for r1 in regions_1.values() {
        for r2 in regions_2.values() {
            // don't intersect regions of different origins
            if r1.origin != r2.origin {
                continue;
            }

            // find the shared substitutions
            let subs_intersect = r1
                .substitutions
                .iter()
                .filter(|sub| r2.substitutions.contains(sub))
                .map(|sub| sub.to_owned())
                .collect::<Vec<_>>();

            // if no shared subs, an intersection is not possible
            if subs_intersect.is_empty() {
                continue;
            }

            // start coordinate is the min sub, end is the max sub
            let start = subs_intersect.iter().min().map(|sub| sub.coord).unwrap();
            let end = subs_intersect.iter().max().map(|sub| sub.coord).unwrap();

            let region = Region {
                start,
                end,
                origin: r1.origin.clone(),
                substitutions: subs_intersect,
            };
            regions_intersect.insert(start, region);
        }
    }

    // Do we need to go back the other way at all?

    Ok(regions_intersect)
}

pub fn identify_breakpoints(
    regions: &BTreeMap<usize, Region>,
) -> Result<Vec<Breakpoint>, Report> {
    let mut breakpoints: Vec<Breakpoint> = Vec::new();
    let mut end_prev: Option<usize> = None;

    for region in regions.values() {
        // a breakpoint is only possible if we already found a prev region
        if let Some(end_prev) = end_prev {
            // breakpoint intervals are non-inclusive of regions
            let breakpoint = Breakpoint {
                start: end_prev + 1,
                end: region.start - 1,
            };
            breakpoints.push(breakpoint);
        }

        end_prev = Some(region.end);
    }

    Ok(breakpoints)
}
