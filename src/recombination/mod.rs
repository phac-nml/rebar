use crate::dataset::{Dataset, SearchResult};
use crate::sequence::Sequence;
use crate::sequence::Substitution;
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

impl Default for Recombination {
    fn default() -> Self {
        Self::new()
    }
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

// ----------------------------------------------------------------------------
// Find Parents

pub struct FindParentsArgs<'a> {

    dataset: &'a Dataset,
    sequence: Sequence,
    best_match: &'a SearchResult,
    max_parents: usize,
    max_iter: usize,
    min_consecutive: usize,
    min_length: usize,
    min_subs: usize,
    
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------


// pub fn find_parents(
//     dataset: &Dataset,
//     sequence: Sequence,
//     best_match: &SearchResult,
//     max_parents: usize,
//     max_iter: usize,
//     min_consecutive: usize,
//     min_length: usize,
//     min_subs: usize,
// ) -> Result<Vec<SearchResult>, Report> {
//     let mut parents = Vec::new();

//     if max_parents == 0 {
//         return Ok(parents);
//     }

//     let mut num_parents = 0;
//     // by default, include all dataset populations
//     let mut include_populations = dataset
//         .populations
//         .keys()
//         .map(|pop| pop.to_owned())
//         .collect::<Vec<_>>();
//     // by default, include all sequence substitutions
//     let mut include_substitutions = sequence.substitutions.clone();
//     // by default, don't exclude any subs or populations
//     let mut exclude_populations = Vec::new();
//     let mut exclude_substitutions: Vec<Substitution> = Vec::new();

//     // --------------------------------------------------------------------
//     // Parent 1
//     // --------------------------------------------------------------------

//     // Option 1. If this is a known recombinant, exclude the recombinant's
//     //   descendants from parent 1 search.

//     debug!("parent_1");
//     let recombinant = best_match.recombinant.clone();
//     let parent_search_result = if let Some(recombinant) = recombinant {
//         let descendants = dataset.phylogeny.get_descendants(&recombinant)?;
//         exclude_populations.extend(descendants);
//         include_populations = include_populations
//             .iter()
//             .filter(|pop| !exclude_populations.contains(&pop))
//             .map(|pop| pop.to_owned())
//             .collect::<Vec<_>>();
//         dataset.search(
//             &sequence,
//             Some(&include_populations),
//             Some(&include_substitutions),
//         )?
//     }
//     // Option 2. Not a known recombinant, just use best_match/consensus as parent 1
//     else {
//         best_match.to_owned()
//     };

//     parents.push(parent_search_result);
//     num_parents += 1;

//     // --------------------------------------------------------------------
//     // Parents 2-MAX
//     // --------------------------------------------------------------------

//     let mut num_iter = 0;

//     loop {
//         // --------------------------------------------------------------------
//         // Loop Break Checks

//         if num_parents >= max_parents {
//             debug!("Maximum number of parents reached ({max_parents}), stopping parent search.");
//             break;
//         }
//         if num_iter >= max_iter {
//             debug!("Maximum number of iterations reached ({num_iter}), stopping parent search.");
//             break;
//         }

//         num_iter += 1;

//         // --------------------------------------------------------------------
//         // Current Mutation Conflicts

//         // first identify all substitutions in all parents
//         let mut parent_substitutions = parents
//             .iter()
//             .map(|p| p.substitutions.clone())
//             .flatten()
//             .collect_vec();
//         parent_substitutions =
//             parent_substitutions.iter().unique().cloned().collect_vec();

//         // identify conflict_alt that are not resolved by another other parent
//         let mut conflict_alt = parents
//             .iter()
//             .map(|p| p.conflict_alt[&p.consensus_population].clone())
//             .flatten()
//             .filter(|sub| !parent_substitutions.contains(&sub))
//             .collect_vec();
//         conflict_alt = conflict_alt.iter().unique().cloned().collect_vec();

//         // conflict_ref is strange, if any parent DOESN'T HAVE the sub, that means it's resolved
//         // not entirely ideal, what about missing and deletions?
//         let mut conflict_ref = parents
//             .iter()
//             .map(|p| p.conflict_ref[&p.consensus_population].clone())
//             .flatten()
//             .collect_vec();
//         conflict_ref = conflict_ref.iter().unique().cloned().collect_vec();

//         let mut conflict_ref_resolved = Vec::new();
//         for sub in &conflict_ref {
//             let is_resolved = false;
//             for parent in &parents {
//                 if !parent.substitutions.contains(sub) {
//                     conflict_ref_resolved.push(sub.to_owned());
//                 }
//             }
//         }

//         conflict_ref = conflict_ref
//             .into_iter()
//             .filter(|sub| !conflict_ref_resolved.contains(sub))
//             .collect_vec();

//         // --------------------------------------------------------------------
//         // Loop Break Checks

//         if conflict_ref.is_empty() {
//             debug!("Sufficient conflict_ref resolution reached, stopping parent search.");
//             break;
//         }
//         if conflict_alt.len() < min_subs {
//             debug!("Sufficient conflict_alt resolution reached, stopping parent search.");
//             break;
//         }

//         debug!("parent_{}: iteration {num_iter}", num_parents + 1);
//         debug!("conflict_ref: {}", conflict_ref.iter().join(", "));
//         debug!("conflict_alt: {}", conflict_alt.iter().join(", "));

//         // --------------------------------------------------------------------
//         // Current Mutation Support

//         let mut parent_support = parents
//             .iter()
//             .map(|p| p.support[&p.consensus_population].clone())
//             .flatten()
//             .collect_vec();
//         parent_support = parent_support.iter().unique().cloned().collect_vec();
//         debug!("sequence_resolved: {}", parent_support.iter().join(", "));

//         let resolved = sequence
//             .substitutions
//             .iter()
//             .filter(|sub| !parent_support.contains(&sub))
//             .cloned()
//             .collect_vec();
//         debug!("sequence_unresolved: {}", resolved.iter().join(", "));

//         // --------------------------------------------------------------------
//         // Filters (include/exclude)

//         // exclude descendants and ancestors(?) of all previous parents
//         for parent in &parents {
//             // descendants to exclude
//             let descendants = dataset
//                 .phylogeny
//                 .get_descendants(&parent.consensus_population)?;
//             let descendants_to_exclude = descendants
//                 .into_iter()
//                 .filter(|pop| !exclude_populations.contains(pop))
//                 .collect::<Vec<_>>();
//             exclude_populations.extend(descendants_to_exclude);

//             // ancestors to exclude
//             let ancestors = dataset
//                 .phylogeny
//                 .get_ancestors(&parent.consensus_population)?;
//             // because of possible recombination, ancestors is a vector of vectors
//             // to reflect multiple parents and paths to the root.
//             // just flatten them for all our purposes here
//             let ancestors_to_add = ancestors
//                 .into_iter()
//                 .flatten()
//                 .filter(|pop| !exclude_populations.contains(pop))
//                 .collect::<Vec<_>>();
//             exclude_populations.extend(ancestors_to_add);
//         }

//         // we want to EXCLUDE populations in our search that
//         //   - have ALL of the conflict_ref
//         let mut conflict_ref_count = BTreeMap::new();
//         for sub in &conflict_ref {
//             if dataset.mutations.contains_key(sub) {
//                 let populations = &dataset.mutations[sub];
//                 for pop in populations {
//                     *conflict_ref_count.entry(pop).or_insert(0) += 1
//                 }
//             }
//         }
//         let populations_to_exclude = conflict_ref_count
//             .into_iter()
//             .filter(|(pop, count)| {
//                 *count == conflict_ref.len()
//                     && !exclude_populations.contains(pop)
//                     && !include_populations.contains(pop)
//             })
//             .map(|(pop, _count)| pop.to_owned())
//             .collect::<Vec<_>>();
//         exclude_populations.extend(populations_to_exclude);

//         // we want to INCLUDE populations in our search that:
//         //   - have AT LEAST min_subs conflict_alt

//         // count up the number of conflict_alt by population
//         let mut conflict_alt_count = BTreeMap::new();
//         for sub in &conflict_alt {
//             if dataset.mutations.contains_key(sub) {
//                 let populations = &dataset.mutations[sub];
//                 for pop in populations {
//                     *conflict_alt_count.entry(pop).or_insert(0) += 1
//                 }
//             }
//         }

//         // reset the include list to just these populations
//         include_populations = conflict_alt_count
//             .into_iter()
//             .filter(|(pop, count)| {
//                 *count >= min_subs && !exclude_populations.contains(pop)
//             })
//             .map(|(pop, _count)| pop.to_owned())
//             .collect::<Vec<_>>();

//         // exclude substitutions that we've already found support for?
//         //include_substitutions = include_substitutions.iter().filter(|sub| )

//         // // reset the include substitutions list
//         // debug!("conflict_ref: {conflict_ref:?}");
//         // include_substitutions = conflict_ref.iter().cloned().cloned().collect::<Vec<_>>();
//         // include_substitutions.extend(conflict_alt);
//         // // include_substitutions = conflict_ref.extend(conflict_alt);

//         // DEBUG
//         //include_populations = vec!(String::from("BA.2.75"));

//         let search_result = dataset.search(
//             &sequence,
//             Some(&include_populations),
//             Some(&include_substitutions),
//         )?;
//         let recombination = detect_recombination(
//             &sequence,
//             &parents,
//             &search_result,
//             min_consecutive,
//             min_length,
//             min_subs,
//         )?;

//         // if the recombination search failed, exclude search_result top_populations from next iteration
//         if recombination.breakpoints.is_empty() {
//             let mut populations_to_exclude = search_result
//                 .top_populations
//                 .iter()
//                 .filter(|pop| !exclude_populations.contains(pop))
//                 .map(|pop| pop.to_owned())
//                 .collect::<Vec<_>>();
//             exclude_populations.append(&mut populations_to_exclude);
//         } else {
//             num_parents += 1;
//             parents.push(search_result);
//         }
//     }

//     Ok(parents)
// }

// pub fn find_breakpoints(
//     sequence: &Sequence,
//     parents: &Vec<dataset::SearchResult>,
//     match_summary: &dataset::SearchResult,
//     min_consecutive: usize,
//     min_length: usize,
//     min_subs: usize,
// ) -> Result<Recombination, Report> {
//     let mut recombination = Recombination::new();

//     // Create a table where rows are coordinates and columns are
//     // coord, parent, Reference, <parents...>, <match_summary> <sequence>

//     let mut table_rows = Vec::new();

//     // Identify which subs are non-bi-allelic, these will wind up being
//     // duplicate rows, which we'll need to reconcile and collapse
//     let mut all_subs = Vec::new();
//     for parent in parents {
//         all_subs.extend(parent.substitutions.to_owned());
//     }
//     all_subs.extend(match_summary.substitutions.to_owned());
//     all_subs.sort();

//     // --------------------------------------------------------------------
//     // Identify Substitution Parental Origins
//     // --------------------------------------------------------------------

//     let all_coords = all_subs.iter().map(|sub| sub.coord).collect::<Vec<_>>();
//     // Deduplicate coords (might be non-bi-allelic)
//     let coords = all_coords.into_iter().unique().collect::<Vec<_>>();
//     let mut privates = Vec::new();

//     for coord in coords {
//         // add coord as first column
//         let mut row = vec![coord.to_string()];

//         // init base origins (could be multiple)
//         let mut origins = Vec::new();

//         // add reference base as second column
//         let ref_base = all_subs
//             .iter()
//             .filter(|sub| sub.coord == coord)
//             .map(|sub| sub.reference)
//             .next()
//             .unwrap();
//         row.push(ref_base.to_string());

//         // Base in Sample
//         let sample_base = sequence
//             .substitutions
//             .iter()
//             .filter(|sub| sub.coord == coord)
//             .map(|sub| sub.alt)
//             .next()
//             .unwrap_or(ref_base);

//         // Add the known parents as next columns
//         let mut parent_bases = Vec::new();

//         for parent in parents {
//             let parent_base = parent
//                 .substitutions
//                 .iter()
//                 .filter(|sub| sub.coord == coord)
//                 .map(|sub| sub.alt)
//                 .next()
//                 .unwrap_or(ref_base);

//             if parent_base == sample_base {
//                 origins.push(parent.consensus_population.clone());
//             }

//             row.push(parent_base.to_string());
//             parent_bases.push(parent_base);
//         }

//         // Add the current parent match that is being evaluated
//         let match_base = match_summary
//             .substitutions
//             .iter()
//             .filter(|sub| sub.coord == coord)
//             .map(|sub| sub.alt)
//             .next()
//             .unwrap_or(ref_base);

//         if match_base == sample_base {
//             origins.push(match_summary.consensus_population.clone());
//         }
//         row.push(match_base.to_string());
//         parent_bases.push(match_base);

//         // Add the sample base
//         row.push(sample_base.to_string());

//         // Is this coord a discriminating site?
//         // Remove subs that are identical between all parents
//         parent_bases = parent_bases.into_iter().unique().collect();

//         // If only 1 unique parent base was found, non-discriminating
//         if parent_bases.len() == 1 {
//             continue;
//         }

//         // If no origins were found, this is private
//         if origins.is_empty() {
//             let private = sequence.substitutions.iter().find(|sub| sub.coord == coord);
//             privates.push(private);
//         } else if origins.len() == 1 {
//             row.push(origins[0].clone());
//             table_rows.push(row);
//         }
//     }

//     // --------------------------------------------------------------------
//     // Debugging Table (pre-filter)
//     // --------------------------------------------------------------------

//     let mut table_builder = Builder::default();

//     // table headers
//     let mut headers = vec!["coord".to_string(), "Reference".to_string()];
//     for parent in parents {
//         headers.push(parent.consensus_population.to_owned());
//     }
//     headers.push(match_summary.consensus_population.to_owned());
//     headers.push(sequence.id.to_owned());
//     headers.push("origin".to_string());

//     table_builder.set_header(headers.clone());

//     // table rows
//     for row in &table_rows {
//         table_builder.push_record(row.to_owned());
//     }

//     let mut table = table_builder.build();

//     // pretty print table for debugging
//     debug!(
//         "Recombination table pre-filter:\n{}",
//         table.with(Style::sharp()).to_string()
//     );

//     // --------------------------------------------------------------------
//     // Group Substitutions into Parental Regions
//     // --------------------------------------------------------------------

//     // First: 5' -> 3', filter separately on min_consecutive then min_length
//     let mut regions_5p = identify_regions(&table_rows)?;
//     regions_5p =
//         filter_regions(&regions_5p, Direction::Forward, min_consecutive, 0_usize)?;
//     regions_5p = filter_regions(&regions_5p, Direction::Forward, 0_usize, min_length)?;
//     debug!(
//         "regions_5p: {}",
//         serde_json::to_string(&regions_5p).unwrap()
//     );

//     // Second: 3' -> 5', filter separately on min_consecutive then min_length
//     let mut regions_3p = identify_regions(&table_rows)?;
//     regions_3p =
//         filter_regions(&regions_3p, Direction::Reverse, min_consecutive, 0_usize)?;
//     regions_3p = filter_regions(&regions_3p, Direction::Reverse, 0_usize, min_length)?;
//     debug!(
//         "regions_3p: {}",
//         serde_json::to_string(&regions_3p).unwrap()
//     );

//     // Take the intersect of the 5' and 3' regions (ie. where they both agree)
//     let mut regions_intersect = intersect_regions(&regions_5p, &regions_3p)?;
//     // During intersection, it's possible that a region from 1 single parent
//     // got broken up into multiple adjacent sections. Put it through the
//     // filter again to collapse it. Direction doesn't matter now
//     regions_intersect = filter_regions(
//         &regions_intersect,
//         Direction::Forward,
//         min_consecutive,
//         min_length,
//     )?;
//     debug!(
//         "regions_intersect: {}",
//         serde_json::to_string(&regions_intersect).unwrap()
//     );

//     // Make sure all the prev_parents + match_summary have at least 1 region
//     let mut region_origins = regions_intersect
//         .values()
//         .map(|region| region.origin.clone());
//     for parent in parents {
//         if !region_origins.contains(&parent.consensus_population) {
//             debug!(
//                 "No recombination detected for parent {}.",
//                 &parent.consensus_population
//             );
//             return Ok(recombination);
//         }
//     }
//     if !region_origins.contains(&match_summary.consensus_population) {
//         debug!(
//             "No recombination detected for parent {}.",
//             &match_summary.consensus_population
//         );
//         return Ok(recombination);
//     }

//     // --------------------------------------------------------------------
//     // Minimum Substitutions Filter
//     // --------------------------------------------------------------------

//     // Check for the minimum number of subs from each parent
//     let mut uniq_subs_count: BTreeMap<String, usize> = BTreeMap::new();
//     // first count up uniq subs by parent
//     for region in regions_intersect.values() {
//         let origin = &region.origin;
//         uniq_subs_count.entry(origin.clone()).or_insert(0);
//         for sub in &region.substitutions {
//             if sub.reference == sub.alt {
//                 continue;
//             }
//             *uniq_subs_count.entry(origin.clone()).or_default() += 1;
//         }
//     }
//     debug!("uniq_subs_count: {uniq_subs_count:?}");
//     // check if any parent failed the filter
//     let mut min_subs_fail = false;
//     for (parent, count) in uniq_subs_count {
//         if count < min_subs {
//             debug!(
//                 "Parent {} subs ({count}) do not meet the min_subs filter ({min_subs}).",
//                 parent
//             );
//             min_subs_fail = true;
//         }
//     }
//     if min_subs_fail {
//         debug!("No recombination detected, min_subs filter was not satisfied by all parents.");
//         return Ok(recombination);
//     }

//     // --------------------------------------------------------------------
//     // Debugging Table (post-filter)
//     // --------------------------------------------------------------------

//     // filter the table for the region subs
//     let mut region_sub_coords = Vec::new();
//     for region in regions_intersect.values() {
//         let sub_coords = region
//             .substitutions
//             .iter()
//             .map(|sub| sub.coord)
//             .collect::<Vec<_>>();
//         region_sub_coords.extend(sub_coords);
//     }

//     let mut table_rows_filter = Vec::new();
//     for row in table_rows {
//         let coord = row[0].parse::<usize>().unwrap();
//         if region_sub_coords.contains(&coord) {
//             table_rows_filter.push(row);
//         }
//     }

//     let mut table_builder = Builder::default();
//     table_builder.set_header(headers.clone());

//     // table rows
//     for row in &table_rows_filter {
//         table_builder.push_record(row.to_owned());
//     }

//     let mut table = table_builder.build();

//     // pretty print table for debugging
//     debug!(
//         "Recombination table post-filter:\n{}",
//         table.with(Style::sharp()).to_string()
//     );

//     // --------------------------------------------------------------------
//     // Breakpoints
//     // --------------------------------------------------------------------

//     let breakpoints = identify_breakpoints(&regions_intersect)?;
//     debug!(
//         "breakpoints: {}",
//         serde_json::to_string(&breakpoints).unwrap()
//     );

//     // --------------------------------------------------------------------
//     // Update
//     // --------------------------------------------------------------------

//     // for the table_rows, combine header and data
//     let mut table = vec![headers];
//     table.extend(table_rows_filter);

//     // update all the attributes
//     recombination.regions = regions_intersect;
//     recombination.breakpoints = breakpoints;
//     recombination.table = table;

//     Ok(recombination)
// }

// pub fn identify_regions(
//     table_rows: &Vec<Vec<String>>,
// ) -> Result<BTreeMap<usize, Region>, Report> {
//     let mut origin_prev: Option<String> = None;
//     let mut regions = BTreeMap::new();
//     let mut start = 0;

//     for row in table_rows {
//         let coord = row[0].parse::<usize>().unwrap();
//         let origin = row[row.len() - 1].to_string();
//         let reference = row[1].chars().next().unwrap();
//         let alt = row[row.len() - 2].chars().next().unwrap();
//         let substitutions = vec![Substitution {
//             coord,
//             reference,
//             alt,
//         }];

//         // start of new region, either first or origin changes
//         if origin_prev.is_none() || origin_prev != Some(origin.clone()) {
//             start = coord;
//             let region = Region {
//                 start,
//                 end: coord,
//                 origin: origin.clone(),
//                 substitutions,
//             };
//             regions.insert(start, region);
//         }
//         // same origin, region continues. update end and subs
//         else {
//             let region = regions.get_mut(&start).unwrap();
//             region.end = coord;
//             region.substitutions.extend(substitutions)
//         }

//         origin_prev = Some(origin);
//     }

//     Ok(regions)
// }

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
