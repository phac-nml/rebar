pub mod parent_search;

use crate::dataset::SearchResult;
use crate::sequence::{Sequence, Substitution};
use color_eyre::eyre::{eyre, Report, Result};
use csv;
use itertools::Itertools;
use log::debug;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::default::Default;
use std::path::Path;
use tabled::builder::Builder;
use tabled::settings::Style;

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
    pub table: Vec<Vec<String>>,
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
            table: Vec::new(),
            genome_length: 0,
        }
    }

    pub fn write_tsv(&self, output_path: &Path) -> Result<(), Report> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(output_path)?;

        for row in &self.table {
            writer.write_record(row)?;
        }

        Ok(())
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

    //recombination.sequence_id = sequence.id.to_string();
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

    let mut table_rows = Vec::new();

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
        // init row with columns:
        // 5 for coord, Reference, Origin, search_result, sequence
        // then variable number of parents
        let num_cols = 5 + parents.len();
        let mut row = vec![String::new(); num_cols];

        // add coord as first column (i=0)
        row[0] = coord.to_string();

        // init base origins (could be multiple)
        // this will be the second columnn (i=1)
        let mut origins = Vec::new();

        // add reference base as third column (i=2)
        let ref_base = all_subs
            .iter()
            .filter(|sub| sub.coord == coord)
            .map(|sub| sub.reference)
            .next()
            .unwrap();
        row[2] = ref_base.to_string();

        // Base in Sample, will use for parent origin comparison
        let sample_base = args
            .sequence
            .substitutions
            .iter()
            .filter(|sub| sub.coord == coord)
            .map(|sub| sub.alt)
            .next()
            .unwrap_or(ref_base);

        // Add the known parents as next columns
        let mut parent_bases = Vec::new();

        for (i, parent) in parents.iter().enumerate() {
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

            row[i + 3] = parent_base.to_string();
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

        // Add search result match_base, second last column
        row[num_cols - 2] = match_base.to_string();
        parent_bases.push(match_base);

        // Add the sample base, last column
        row[num_cols - 1] = sample_base.to_string();

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
            // origins is second columnn (i=1)
            row[1] = origins.iter().join(",");
            table_rows.push(row);
        }
    }

    // --------------------------------------------------------------------
    // Group Substitutions into Parental Regions
    // --------------------------------------------------------------------

    // First: 5' -> 3', filter separately on min_consecutive then min_length
    let mut regions_5p = identify_regions(&table_rows)?;
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
    let mut regions_3p = identify_regions(&table_rows)?;
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
                parent, args.min_subs,
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
    let mut headers = vec!["coord", "origin", "Reference"]
        .into_iter()
        .map(String::from)
        .collect_vec();

    for parent in parents {
        headers.push(parent.consensus_population.to_owned());
    }
    headers.push(search_result.consensus_population.to_owned());
    headers.push(args.sequence.id.to_owned());

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
    recombination.parents = region_origins;
    recombination.regions = regions_intersect;
    recombination.breakpoints = breakpoints;
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
    table_rows: &Vec<Vec<String>>,
) -> Result<BTreeMap<usize, Region>, Report> {
    let mut origin_prev: Option<String> = None;
    let mut regions = BTreeMap::new();
    let mut start = 0;

    for row in table_rows {
        let coord = row[0].parse::<usize>().unwrap();
        let origin = row[1].to_string();
        let reference = row[2].chars().next().unwrap();
        let alt = row[row.len() - 1].chars().next().unwrap();
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
) -> Result<Vec<Vec<String>>, Report> {
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
    let mut combine_table: Vec<Vec<String>> = Vec::new();

    // Mandatory headers (coord and Reference)
    let mut headers = vec!["coord", "origin", "Reference"];
    // Dynamic headers (parents and sequence IDs)
    for parent in parents.iter() {
        headers.push(parent)
    }
    for sequence_id in sequence_ids.iter() {
        headers.push(sequence_id)
    }
    // convert to String, &str won't work here, since we're going to create
    // table row values within a for loop scope later
    let headers = headers.into_iter().map(String::from).collect_vec();
    combine_table.push(headers.clone());

    // identify all coords in all samples, convert from string to numeric
    // so they can be sorted nicely
    let mut coords = recombinations
        .iter()
        .flat_map(|rec| rec.table.iter().map(|row| &row[0]).collect_vec())
        .unique()
        .filter(|c| c.as_str() != "coord")
        .map(|c| c.parse::<usize>().unwrap())
        .collect_vec();
    coords.sort();

    // combine tables, row by row (coord by coord)
    for coord in &coords {
        // init row with empty strings for all columns
        let mut row = vec![String::new(); headers.len()];
        // get reference base directly from sequence
        let ref_base = reference.seq[coord - 1].to_string();
        // Reference is always the third col (i=2)
        row[2] = ref_base.to_string();

        // iterate through recombinants, identifying ref, parents, seq bases
        for (rec_i, recombination) in recombinations.iter().enumerate() {
            // get sequence base directly from sequence
            let rec_base = recombination.sequence.seq[coord - 1].to_string();
            // what is the col position of this sample?
            // first 3 are coord, origin, ref, then parents, then recs
            let col_i = 3 + parents.len() + rec_i;
            row[col_i] = rec_base;

            // store the recombinant barcode table
            // skip the first row (header)
            let rec_table = &recombination.table.iter().skip(1).collect_vec();

            // check which coords appeared in this sample
            let rec_coords = rec_table
                .iter()
                .map(|row| row[0].parse::<usize>().unwrap())
                .collect_vec();

            // get table index of coord
            let row_i = rec_coords.iter().position(|c| c == coord);

            // if this sample has the coord
            if let Some(row_i) = row_i {
                // if it's the first sample encountered add the ref and parent bases
                if row[0] == String::new() {
                    // Add the coord to the table row, coord is always first col (i=0)
                    row[0] = coord.to_string();

                    // Origin is second col (i=1)
                    let origin = &rec_table[row_i][1];
                    row[1] = origin.to_string();

                    // Add parents, i + 3, after coord (0), Origin (1), Reference (2)
                    for (i, _parent) in parents.iter().enumerate() {
                        let parent_base = &rec_table[row_i][i + 3];
                        row[i + 3] = parent_base.to_string();
                    }
                }
            }

            // let rec_base = recombinant.sequence.seq

            // // convert deletions to just coords for checking
            // let deletion_coords = recombination
            //     .sequence
            //     .deletions
            //     .iter()
            //     .map(|del| del.coord)
            //     .collect_vec();

            // // convert subs to just coords for checking
            // let substitution_coords = recombination
            //     .sequence
            //     .substitutions
            //     .iter()
            //     .map(|sub| sub.coord)
            //     .collect_vec();

            // // missing
            // if recombination.sequence.missing.contains(coord) {
            //     rec_base = "N";
            // }
            // // substitution
            // else if substitutions_coords.contains(coord) {
            //     rec_base = recombination.sequence.substitutions
            // }
            // // deletion
            // else if deletion_coords.contains(coord) {
            //     rec_base = "-";
            // }

            // what is the col position of this sample?
            // first 3 are coord, origin, ref, then parents, then recs
            //let col_i = 3 + parents.len() + rec_i;
            //row[col_i] = rec_base.to_string();
        }

        // Add processed row to table
        combine_table.push(row);
    }

    Ok(combine_table)
}
