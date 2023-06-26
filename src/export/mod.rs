use crate::dataset::{Dataset, SearchResult};
use crate::recombination::{combine_barcode_tables, Recombination};
use color_eyre::eyre::{eyre, Report, Result};
use csv;
use itertools::Itertools;
use log::debug;
//use log::debug, info};
//use std::collections::BTreeMap;
use std::path::Path;

// ----------------------------------------------------------------------------
// LineList

#[derive(Debug)]
pub struct LineList {
    strain: Vec<String>,
    population: Vec<String>,
    recombinant: Vec<String>,
    parents: Vec<String>,
    breakpoints: Vec<String>,
    unique_key: Vec<String>,
    regions: Vec<String>,
    private: Vec<String>,
    genome_length: Vec<String>,
    dataset_name: Vec<String>,
    dataset_tag: Vec<String>,
}

impl Default for LineList {
    fn default() -> Self {
        Self::new()
    }
}

impl LineList {
    pub fn new() -> Self {
        LineList {
            strain: Vec::new(),
            population: Vec::new(),
            recombinant: Vec::new(),
            parents: Vec::new(),
            breakpoints: Vec::new(),
            unique_key: Vec::new(),
            regions: Vec::new(),
            private: Vec::new(),
            genome_length: Vec::new(),
            dataset_name: Vec::new(),
            dataset_tag: Vec::new(),
        }
    }

    // Get LineList attribute by str.
    pub fn get(&self, attribute: &str) -> Result<&Vec<String>, Report> {
        let result = match attribute {
            "strain" => &self.strain,
            "population" => &self.population,
            "recombinant" => &self.recombinant,
            "parents" => &self.parents,
            "breakpoints" => &self.breakpoints,
            "unique_key" => &self.unique_key,
            "regions" => &self.regions,
            "private" => &self.private,
            "genome_length" => &self.genome_length,
            "dataset_name" => &self.dataset_name,
            "dataset_tag" => &self.dataset_tag,
            _ => return Err(eyre!("Attribute {attribute} is not in LineList.")),
        };

        Ok(result)
    }

    pub fn headers() -> Vec<String> {
        vec![
            "strain",
            "population",
            "recombinant",
            "parents",
            "breakpoints",
            "unique_key",
            "regions",
            "private",
            "genome_length",
            "dataset_name",
            "dataset_tag",
        ]
        .into_iter()
        .map(|s| s.to_string())
        .collect_vec()
    }

    pub fn write_tsv(&self, output_path: &Path) -> Result<(), Report> {
        // initialize writer
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(output_path)?;

        // write headers
        let headers = LineList::headers();
        writer.write_record(&headers)?;

        let num_rows = self.strain.len();

        // if the linelist was empty, just write headers
        if num_rows == 0 {
            return Ok(());
        }

        for row_i in 0..num_rows {
            let mut row = vec![""; headers.len()];
            for (col_i, header) in headers.iter().enumerate() {
                let vals = self.get(header)?;
                row[col_i] = &vals[row_i];
            }
            writer.write_record(row)?;
        }

        Ok(())
    }

    pub fn create(
        recombinations: &Vec<Recombination>,
        best_matches: &Vec<SearchResult>,
        dataset: &Dataset,
    ) -> Result<LineList, Report> {
        // init linelist
        let mut linelist = LineList::new();

        // check for same length
        if recombinations.len() != best_matches.len() {
            return Err(eyre!(
                "recombinations and best_matches are different lengths."
            ));
        }

        // iterate in parallel, checking for same sequence id
        for it in recombinations.iter().zip(best_matches.iter()) {
            let (recombination, best_match) = it;

            // check that they're in the correct order
            if recombination.sequence.id != best_match.sequence_id {
                return Err(eyre!(
                    "Recombination ID {} is not the same as Best Match ID: {}",
                    recombination.sequence.id,
                    best_match.sequence_id,
                ));
            }

            // strain
            let strain = recombination.sequence.id.to_string();
            linelist.strain.push(strain);

            // population
            let population = best_match.consensus_population.to_string();
            linelist.population.push(population);

            // recombinant
            let recombinant = &best_match.recombinant;
            match recombinant {
                Some(recombinant) => linelist.recombinant.push(recombinant.to_string()),
                None => linelist.recombinant.push(String::new()),
            }

            // parents
            let parents = recombination.parents.join(",").to_string();
            linelist.parents.push(parents);

            // breakpoints
            let breakpoints = recombination.breakpoints.iter().join(",").to_string();
            linelist.breakpoints.push(breakpoints);

            // unique_key
            let unique_key = recombination.unique_key.to_string();
            linelist.unique_key.push(unique_key);

            // regions
            let regions = recombination.regions.values().join(",").to_string();
            linelist.regions.push(regions);

            // private mutations
            let private = best_match.private.iter().join(",").to_string();
            linelist.private.push(private);

            // genome_length
            let genome_length = recombination.genome_length.to_string();
            linelist.genome_length.push(genome_length);

            // dataset name
            linelist.dataset_name.push(dataset.name.to_string());

            // dataset tag
            linelist.dataset_tag.push(dataset.tag.to_string());
        }

        Ok(linelist)
    }
}

pub fn write_barcodes(
    _output_dir: &Path,
    linelist: &LineList,
    recombinations: &[Recombination],
) -> Result<(), Report> {
    let num_rows = linelist.strain.len();
    let unique_keys = linelist.unique_key.iter().unique().collect_vec();
    debug!("{unique_keys:?}");

    for unique_key in unique_keys {
        let mut strains = Vec::new();
        // identify strains part belong to this key
        for i in 0..num_rows {
            if &linelist.unique_key[i] == unique_key {
                strains.push(&linelist.strain[i]);
            }
        }

        // identify recombinations part of this key
        let unique_rec = recombinations
            .iter()
            .filter(|rec| strains.contains(&&rec.sequence.id))
            .cloned()
            .collect_vec();

        // combine recombination barcode tables
        combine_barcode_tables(&unique_rec)?;
    }

    Ok(())
}
