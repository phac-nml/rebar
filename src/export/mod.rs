use crate::dataset::{Dataset, SearchResult};
use crate::recombination::Recombination;
use color_eyre::eyre::{eyre, Report, Result};
use csv;
use itertools::Itertools;
use log::{debug, info};
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Debug)]
pub struct LineList {
    strain: Vec<String>,
    population: Vec<String>,
    recombinant: Vec<String>,
    parents: Vec<String>,
    breakpoints: Vec<String>,
    regions: Vec<String>,
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
            regions: Vec::new(),
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
            "regions" => &self.regions,
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
            "regions",
            "genome_length",
            "dataset_name",
            "dataset_tag",
        ]
        .into_iter()
        .map(|s| s.to_string())
        .collect_vec()
    }

    /// Collect strains that belong to the same unique recombinant.
    ///
    /// A unique recombinant is defined by having the same:
    ///   - Recombinant population
    ///   - Parents
    ///   - Breakpoints
    pub fn collect_recombinants(&self) -> Result<(), Report> {
        let mut uniq_recombinants: BTreeMap<String, Vec<String>> = BTreeMap::new();

        // Collect recombinants into map where keys are <recombinant_<parents>_<breakpoints>
        // and values are Vec of strain names.

        let num_rows = self.strain.len();
        for i in 0..num_rows {
            let strain = &self.strain[i];
            let recombinant = &self.recombinant[i];
            let parents = &self.parents[i];
            let breakpoints = &self.breakpoints[i];

            let uniq_key = format!("{}_{}_{}", recombinant, parents, breakpoints);

            uniq_recombinants
                .entry(uniq_key.to_string())
                .or_insert(Vec::new())
                .push(strain.to_string());
        }

        // Iterate through the collection
        for (uniq_key, strains) in uniq_recombinants {
            info!("{uniq_key}: {strains:?}");

            // combine barcodes
        }

        Ok(())
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

    pub fn from_recombinations(
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

        // iterate in parallel, checking for same sequence_id
        for it in recombinations.iter().zip(best_matches.iter()) {
            let (recombination, best_match) = it;

            // check that they're in the correct order
            if recombination.sequence_id != best_match.sequence_id {
                return Err(eyre!(
                    "Recombination ID {} is not the same as Best Match ID: {}",
                    recombination.sequence_id,
                    best_match.sequence_id,
                ));
            }

            // strain
            let strain = recombination.sequence_id.to_string();
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

            // regions
            let regions = recombination.regions.values().join(",").to_string();
            linelist.regions.push(regions);

            // genome_length
            let genome_length = recombination.genome_length.to_string();
            linelist.genome_length.push(genome_length);

            // dataset name
            linelist.dataset_name.push(dataset.name.to_string());

            // dataset tag
            linelist.dataset_tag.push(dataset.tag.to_string());
        }

        debug!("{linelist:?}");

        Ok(linelist)
    }
}

// pub fn write_barcodes(
//     output_dir: &Path,
//     sequences: &Vec<Sequence>,
//     recombinations: &BTreeMap<String, Recombination>,
// ) -> Result<(), Report> {
//     Ok(())
// }
