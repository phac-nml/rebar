pub mod parent_search;

use crate::dataset::SearchResult;
use crate::sequence::{Sequence, Substitution};
use crate::utils;
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;
use log::debug;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

// ----------------------------------------------------------------------------
// Structs
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Breakpoint

/// Recombination breakpoint intervals (left and right inclusive)
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Breakpoint {
    pub start: usize,
    pub end: usize,
}

impl std::fmt::Display for Breakpoint {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}-{}", self.start, self.end)
    }
}

// ----------------------------------------------------------------------------
// Direction

/// Genomic reading direction as forward (5' -> 3') or reverse (3' -> 5')
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

impl std::fmt::Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}-{}|{}", self.start, self.end, self.origin)
    }
}

// ----------------------------------------------------------------------------
// Recombination

#[derive(Clone, Debug, Serialize)]
pub struct Recombination<'seq> {
    //pub sequence_id: String,
    pub sequence: &'seq Sequence,
    pub unique_key: String,
    pub recombinant: Option<String>,
    pub parents: Vec<String>,
    pub breakpoints: Vec<Breakpoint>,
    pub regions: BTreeMap<usize, Region>,
    pub genome_length: usize,
    #[serde(skip_serializing)]
    pub table: utils::table::Table,
}

impl<'seq> Recombination<'seq> {
    pub fn new(sequence: &'seq Sequence) -> Self {
        Recombination {
            sequence,
            unique_key: String::new(),
            recombinant: None,
            parents: Vec::new(),
            breakpoints: Vec::new(),
            regions: BTreeMap::new(),
            table: utils::table::Table::new(),
            genome_length: 0,
        }
    }
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

pub fn detect_recombination<'seq>(
    parents: &Vec<SearchResult>,
    search_result: Option<&SearchResult>,
    args: &parent_search::Args<'_, 'seq, '_>,
) -> Result<Recombination<'seq>, Report> {
    let mut recombination = Recombination::new(args.sequence);

    recombination.sequence = args.sequence;
    recombination.genome_length = args.sequence.genome_length;

    // if no search_result was provided, just use first parent
    let search_result_default = &parents[0];
    let search_result = match search_result {
        Some(summary) => summary,
        None => search_result_default,
    };

    // Create a table where rows are coordinates and columns are
    // coord, parent, Reference, <parents...>, <search_result> <sequence>

    let mut table = utils::table::Table::new();
    table.headers = vec!["coord", "origin", "Reference"]
        .into_iter()
        .map(String::from)
        .collect_vec();
    for parent in parents {
        table.headers.push(parent.consensus_population.to_string());
    }
    table
        .headers
        .push(search_result.consensus_population.to_string());
    table.headers.push(args.sequence.id.to_string());
    table.rows = Vec::new();

    // get the column position of each header, trying to avoid hard-coding
    let coord_col_i = table.header_position("coord")?;
    let origin_col_i = table.header_position("origin")?;
    let ref_col_i = table.header_position("Reference")?;
    let search_col_i = table.header_position(&search_result.consensus_population)?;
    let seq_col_i = table.header_position(&args.sequence.id)?;

    // Identify all subs found in parents and search_result
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
        // init row with columns
        let mut row = vec![String::new(); table.headers.len()];

        // add coord to row
        row[coord_col_i] = coord.to_string();

        // get Reference base directly from sequence
        let ref_base = &args.dataset.reference.seq[coord - 1];
        row[ref_col_i] = ref_base.to_string();

        // get Sample base directly form sequence
        let seq_base = &args.sequence.seq[coord - 1];
        row[seq_col_i] = seq_base.to_string();

        // init base origins (could be multiple)
        let mut origins = Vec::new();

        // Keep track of parent bases diversity, if they're all the same, this
        // isn't an informative/differentiating site
        let mut parent_bases = Vec::new();

        for parent in parents {
            let parent_base = parent
                .substitutions
                .iter()
                .filter(|sub| sub.coord == coord)
                .map(|sub| sub.alt)
                .next()
                .unwrap_or(*ref_base);

            if parent_base == *seq_base {
                origins.push(parent.consensus_population.clone());
            }

            let parent_col_i = table.header_position(&parent.consensus_population)?;
            row[parent_col_i] = parent_base.to_string();
            parent_bases.push(parent_base);
        }

        // Add the search result that is currently being evaluated
        let search_result_base = search_result
            .substitutions
            .iter()
            .filter(|sub| sub.coord == coord)
            .map(|sub| sub.alt)
            .next()
            .unwrap_or(*ref_base);

        if search_result_base == *seq_base {
            origins.push(search_result.consensus_population.clone());
        }

        // Add search result to the row
        row[search_col_i] = search_result_base.to_string();
        parent_bases.push(search_result_base);

        // Is this coord a discriminating site?
        // Remove subs that are identical between all parents
        parent_bases = parent_bases.into_iter().unique().collect();

        // If only 1 unique parent base was found, non-discriminating
        if parent_bases.len() == 1 {
            continue;
        }

        // If no origins were found, this is private
        if origins.is_empty() {
            let private = args
                .sequence
                .substitutions
                .iter()
                .find(|sub| sub.coord == coord);
            privates.push(private);
        } else if origins.len() == 1 {
            // add origins to row
            row[origin_col_i] = origins.iter().join(",");
            table.rows.push(row);
        }
    }

    // Debugging Table
    debug!("Recombination table:\n\n{}", table.to_markdown()?);

    // --------------------------------------------------------------------
    // Group Substitutions into Parental Regions
    // --------------------------------------------------------------------

    // First: 5' -> 3', filter separately on min_consecutive then min_length
    let mut regions_5p = identify_regions(&table)?;
    regions_5p = filter_regions(
        &regions_5p,
        Direction::Forward,
        args.min_consecutive,
        0_usize,
    )?;
    regions_5p =
        filter_regions(&regions_5p, Direction::Forward, 0_usize, args.min_length)?;
    debug!(
        "regions_5p: {}",
        serde_json::to_string(&regions_5p).unwrap()
    );

    // Second: 3' -> 5', filter separately on min_consecutive then min_length
    let mut regions_3p = identify_regions(&table)?;
    regions_3p = filter_regions(
        &regions_3p,
        Direction::Reverse,
        args.min_consecutive,
        0_usize,
    )?;
    regions_3p =
        filter_regions(&regions_3p, Direction::Reverse, 0_usize, args.min_length)?;
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
        .unique()
        .collect_vec();

    for parent in parents {
        if !region_origins.contains(&parent.consensus_population) {
            return Err(eyre!(
                "No recombination detected for parent {}.",
                &parent.consensus_population
            ));
        }
    }
    if !region_origins.contains(&search_result.consensus_population) {
        return Err(eyre!(
            "No recombination detected for parent {}.",
            &search_result.consensus_population
        ));
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
                parent, args.min_subs,
            );
            min_subs_fail = true;
        }
    }
    if min_subs_fail {
        return Err(eyre!("No recombination detected, min_subs filter was not satisfied by all parents."));
    }

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

    // update all the attributes
    recombination.parents = region_origins;
    recombination.regions = regions_intersect;
    recombination.breakpoints = breakpoints;
    // pre-filter or post-filter table?
    recombination.table = table;
    recombination.recombinant = args.best_match.recombinant.clone();

    // Grab name of recombinant for unique key
    let best_match_recombinant = &args.best_match.recombinant;
    let rec_key = match best_match_recombinant {
        Some(recombinant) => recombinant,
        None => "novel",
    };

    recombination.unique_key = format!(
        "{}_{}_{}",
        rec_key,
        &recombination.parents.iter().join("_"),
        &recombination.breakpoints.iter().join("_"),
    );

    Ok(recombination)
}

pub fn identify_regions(
    table: &utils::table::Table,
) -> Result<BTreeMap<usize, Region>, Report> {
    let mut origin_prev: Option<String> = None;
    let mut regions = BTreeMap::new();
    let mut start = 0;

    let coord_col_i = table.header_position("coord")?;
    let origin_col_i = table.header_position("origin")?;
    let ref_col_i = table.header_position("Reference")?;
    // sequence for alt is last column
    let seq_col_i = table.headers.len() - 1;

    for row in table.rows.iter() {
        let coord = row[coord_col_i].parse::<usize>().unwrap();
        let origin = row[origin_col_i].to_string();
        let reference = row[ref_col_i].chars().next().unwrap();
        let alt = row[seq_col_i].chars().next().unwrap();
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

/// Identify breakpoint intervals in recombination regions.
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

/// Combine recombination tables.
pub fn combine_tables(
    recombinations: &[Recombination],
    reference: &Sequence,
) -> Result<utils::table::Table, Report> {
    // ------------------------------------------------------------------------
    // Input Checking

    // was recombination detected?
    let parents = recombinations.iter().map(|rec| &rec.parents).collect_vec();
    if parents.is_empty() {
        return Err(eyre!("Cannot combine tables, no recombination detected."));
    }

    // are all parents the same?
    for parent in &parents {
        if parent != &parents[0] {
            return Err(eyre!(
                "Cannot combine tables, not all parents are the same."
            ));
        }
    }

    // identify sequence IDs to combine
    let sequence_ids = recombinations
        .iter()
        .map(|rec| &rec.sequence.id)
        .collect_vec();

    // identify parents to combine, just use first, since we verified all same
    let parents = recombinations
        .iter()
        .map(|rec| &rec.parents)
        .next()
        .unwrap();

    // ------------------------------------------------------------------------
    // Construct Table Headers

    // final result to mutate and return
    let mut combine_table = utils::table::Table::new();

    // Mandatory headers
    // convert to String, &str won't work here, since we're going to create
    // table row values within a for loop scope later
    combine_table.headers = vec!["coord", "origin", "Reference"]
        .into_iter()
        .map(String::from)
        .collect_vec();

    // Dynamic headers (parents and sequence IDs)
    for parent in parents {
        combine_table.headers.push(parent.clone())
    }
    for sequence_id in sequence_ids {
        combine_table.headers.push(sequence_id.clone())
    }

    // get output col idx of mandatory columns
    let coord_output_i = combine_table.header_position("coord")?;
    let origin_output_i = combine_table.header_position("origin")?;
    let ref_output_i = combine_table.header_position("Reference")?;

    // ------------------------------------------------------------------------
    // Combine Coordinate Bases

    // identify all coords in all samples, convert from string to numeric
    // so they can be sorted nicely
    let mut coords = recombinations
        .iter()
        .flat_map(|rec| rec.table.rows.iter().map(|row| &row[0]).collect_vec())
        .unique()
        .filter(|c| c.as_str() != "coord")
        .map(|c| c.parse::<usize>().unwrap())
        .collect_vec();
    coords.sort();

    // combine tables, row by row (coord by coord)
    for coord in &coords {
        // init row with empty strings for all columns
        let mut row = vec![String::new(); combine_table.headers.len()];
        // get reference base directly from sequence
        let ref_base = reference.seq[coord - 1].to_string();
        row[ref_output_i] = ref_base.to_string();

        // iterate through recombinants, identifying ref, parents, seq bases
        for recombination in recombinations {
            // get sequence base directly from sequence
            let rec_base = recombination.sequence.seq[coord - 1].to_string();
            let rec_output_i =
                combine_table.header_position(&recombination.sequence.id)?;
            row[rec_output_i] = rec_base;

            // check which coords appeared in this sample
            let coord_input_i = recombination.table.header_position("coord")?;
            let rec_coords = recombination
                .table
                .rows
                .iter()
                .map(|row| row[coord_input_i].parse::<usize>().unwrap())
                .collect_vec();

            // get table index of these coords
            let row_input_i = rec_coords.iter().position(|c| c == coord);

            // get table index of origin
            let origin_input_i = recombination.table.header_position("origin")?;

            // if this sample has the coord
            if let Some(row_input_i) = row_input_i {
                // if it's the first sample encountered add the coord, origin,
                // Reference base and parent bases
                if row[coord_output_i] == String::new() {
                    // Add the coord to the output table
                    row[coord_output_i] = coord.to_string();

                    // add the origin to the output table
                    let origin = &recombination.table.rows[row_input_i][origin_input_i];
                    row[origin_output_i] = origin.to_string();

                    // Add parents bases to the output table
                    for parent in parents {
                        let parent_input_i =
                            recombination.table.header_position(parent)?;
                        let parent_output_i = combine_table.header_position(parent)?;
                        let parent_base =
                            &recombination.table.rows[row_input_i][parent_input_i];
                        row[parent_output_i] = parent_base.to_string();
                    }
                }
            }
        }

        // Add processed row to table
        combine_table.rows.push(row);
    }

    Ok(combine_table)
}
