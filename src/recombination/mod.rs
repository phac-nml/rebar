use crate::query::match_summary::MatchSummary;
use crate::sequence::Sequence;
use crate::sequence::Substitution;
use color_eyre::eyre::{Report, Result};
use itertools::Itertools;
use log::debug;
use std::collections::BTreeMap;
use std::default::Default;
use tabled::builder::Builder;
use tabled::settings::Style;

#[derive(Debug)]
pub struct Breakpoint {
    start: usize,
    end: usize,
}

#[derive(Debug)]
#[allow(dead_code)]
pub struct Region {
    start: usize,
    end: usize,
    origin: String,
    substitutions: Vec<Substitution>,
}

#[derive(Debug)]
pub struct Recombination {
    breakpoints: Vec<Breakpoint>,
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
        }
    }

    pub fn detect(
        &mut self,
        sequence: &Sequence,
        parents: &Vec<MatchSummary>,
        match_summary: &MatchSummary,
    ) -> Result<(), Report> {
        // Create a table where rows are coordinates and columns are
        // coord, parent, Reference, <parents...>, <match_summary> <sequence>

        let mut table_rows = Vec::new();

        // Identify which subs are non-bi-allelic, these will wind up being
        // duplicate rows, which we'll need to reconcile and collapse
        let mut all_subs = Vec::new();
        for parent in parents {
            all_subs.extend(parent.substitutions.to_owned());
        }

        all_subs.sort();
        let all_coords = all_subs.iter().map(|sub| sub.coord).collect::<Vec<_>>();

        // --------------------------------------------------------------------
        // Identify Substitution Parental Origins
        // --------------------------------------------------------------------

        let mut privates = Vec::new();

        let coords = all_coords.into_iter().unique().collect::<Vec<_>>();

        for coord in coords {
            // add coord as first column
            let mut row = vec![coord.to_string()];

            // init sub origins (could be multiple)
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
            let match_base = match_summary
                .substitutions
                .iter()
                .filter(|sub| sub.coord == coord)
                .map(|sub| sub.alt)
                .next()
                .unwrap_or(ref_base);

            if match_base == sample_base {
                origins.push(match_summary.consensus_population.clone());
            }
            row.push(match_base.to_string());
            parent_bases.push(match_base);

            // Add the sample base
            row.push(sample_base.to_string());

            // Is this coord a discriminating site?
            // Remove subs that are identical between all parents
            parent_bases = parent_bases.into_iter().unique().collect();

            // If no origins were found, this is private
            if origins.is_empty() {
                let private =
                    sequence.substitutions.iter().find(|sub| sub.coord == coord);
                privates.push(private);
                continue;
            }
            // If only 1 unique parent base was found, non-discriminating
            else if parent_bases.len() == 1 {
                continue;
            }

            // sub origins (csv list of populations)
            let origins = origins.iter().join(", ");
            row.push(origins);

            table_rows.push(row);
        }

        // --------------------------------------------------------------------
        // Table
        // --------------------------------------------------------------------

        let mut table_builder = Builder::default();

        // column headers
        let mut headers = vec!["coord".to_string(), "Reference".to_string()];
        for parent in parents {
            headers.push(parent.consensus_population.to_owned());
        }
        headers.push(match_summary.consensus_population.to_owned());
        headers.push(sequence.id.to_owned());
        headers.push("origin".to_string());

        table_builder.set_header(headers);

        for row in &table_rows {
            table_builder.push_record(row.to_owned());
        }

        let mut table = table_builder.build();
        debug!(
            "{}",
            table
                .with(Style::sharp())
                .to_string()
                .replace('\n', format!("\n{}", " ".repeat(46)).as_str())
        );

        // --------------------------------------------------------------------
        // Group Substitutions Into Parental Regions
        // --------------------------------------------------------------------

        // First: 5' -> 3'
        let regions_5p = identify_regions(&table_rows)?;
        // Filter separately on min_consecutive then min_length
        debug!("{regions_5p:?}");
        // Second: 5' -> 3'

        // --------------------------------------------------------------------
        // Breakpoints
        // --------------------------------------------------------------------

        self.breakpoints = vec![Breakpoint {
            start: 123,
            end: 456,
        }];
        debug!("{} {}", self.breakpoints[0].start, self.breakpoints[0].end);

        // create the summary table
        Ok(())
    }
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
        let reference = row[2].chars().next().unwrap();
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
