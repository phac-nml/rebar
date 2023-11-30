use crate::cli;
use crate::dataset;
use crate::recombination;

use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use itertools::Itertools;
use log::{debug, info};
use rand::Rng;
use std::fs::{create_dir_all, File};
use std::io::Write;

/// Simulate recombination.
pub fn simulate(args: &cli::simulate::Args) -> Result<(), Report> {
    // create output directory if it doesn't exist
    if !args.output_dir.exists() {
        info!("Creating output directory: {:?}", args.output_dir);
        create_dir_all(&args.output_dir)?;
    }

    // Load dataset, disable masking
    info!("Loading dataset: {:?}", &args.dataset_dir);
    let mask = vec![0, 0];
    let dataset = dataset::load::dataset(&args.dataset_dir, &mask)?;
    let genome_length = dataset.reference.genome_length;

    // Check to make sure all parents are in dataset
    let parents = args.parents.clone();
    for parent in &parents {
        if !dataset.populations.contains_key(parent.as_str()) {
            return Err(eyre!(
                "Parent {parent} is not the dataset populations fasta."
            ));
        }
    }

    // ------------------------------------------------------------------------
    // Breakpoints

    let mut breakpoints = if let Some(breakpoints) = &args.breakpoints {
        info!("Using manual breakpoints: {breakpoints:?}");
        breakpoints.clone()
    } else {
        let mut breakpoints = Vec::new();
        let mut num_breakpoints_remaining = parents.len() - 1;
        let mut start = 1;
        let mut rng = rand::thread_rng();
        while num_breakpoints_remaining > 0 {
            // save some coordinates for future breakpoints
            let end = genome_length - num_breakpoints_remaining;
            let coord = rng.gen_range(start..end);
            breakpoints.push(coord);
            start = coord + 1;
            num_breakpoints_remaining -= 1;
        }
        info!("Using random breakpoints: {breakpoints:?}");
        breakpoints
    };

    let unique_key = format!(
        "simulate_{}_{}",
        &parents.iter().join("_"),
        &breakpoints.iter().map(|start| format!("{start}-{}", start + 1)).join("_"),
    );
    info!("Unique Key: {unique_key:?}");

    // ------------------------------------------------------------------------
    // Regions

    breakpoints.push(genome_length);
    let mut regions = Vec::new();
    let mut start = 1;

    for (origin, end) in parents.into_iter().zip(breakpoints.into_iter()) {
        let region = recombination::Region {
            start,
            end,
            origin,
            substitutions: Vec::new(),
        };
        regions.push(region);
        start = end + 1;
    }
    debug!("Regions: {regions:?}");

    // ------------------------------------------------------------------------
    // Sequences

    let sequence: String = regions
        .iter()
        .map(|region| {
            let sequence = dataset.populations.get(&region.origin).unwrap_or_else(|| {
                panic!(
                    "Failed to find region origin {} in dataset populations.",
                    &region.origin
                )
            });
            // Reminder, -1 to coordinates since they are 1-based
            sequence.seq[region.start - 1..=region.end - 1].iter().collect::<String>()
        })
        .collect();

    let output_path = args.output_dir.join(format!("{unique_key}.fasta"));
    info!("Exporting fasta: {output_path:?}");
    let mut output_file = File::create(&output_path)
        .wrap_err_with(|| format!("Unable to create file: {output_path:?}"))?;
    let lines = format!(">{unique_key}\n{sequence}");
    output_file
        .write_all(lines.as_bytes())
        .wrap_err_with(|| format!("Unable to write file: {output_path:?}"))?;

    Ok(())
}
