pub mod cli;
pub mod dataset;
pub mod phylogeny;
pub mod query;
pub mod recombination;
pub mod sequence;
pub mod traits;
pub mod utils;

use bio::io::fasta;
use crate::sequence::Sequence;
use cli::RebarRunArgs;
use color_eyre::eyre::{Result, Report, WrapErr};
use log::{debug, info, warn};

/// Run rebar on input alignment or dataset population
pub fn run(args: RebarRunArgs) -> Result<(), Report> {
    println!("{args:?}");

    // Collect files in dataset_dir into a dataset object
    // This mainly includes parent populations sequences
    //   and optionally a phylogenetic representation.
    let dataset = dataset::load(&args.dataset_dir, args.mask)?;  

    // init a containter to hold query sequences, either dataset 
    // populations or sequences from an input alignment
    let mut sequences = Vec::new();

    // ------------------------------------------------------------------------
    // input population was specified

    let populations = args.input.populations;
    if let Some( populations ) = populations {

        // split population into csv (,) parts
        let populations = populations.split(",").collect::<Vec<_>>();

        let mut search_populations = Vec::new();
    
        // if population ends in '*', use phylogenetically aware mode
        // to search for it and all descendants
        for population in populations {

            if population.ends_with("*"){
                // remove last char (wildcard) with pop
                let parent = population[0..population.len() - 1].to_string();
                // get descendants of all parents
                let descendants = dataset.phylogeny.get_descendants(&parent)?;
                // add parent population and descendants
                search_populations.push(parent);
                search_populations.extend(descendants);
            } else {
                search_populations.push(population.to_string());
            }
        }

        for population in search_populations {        
            // check if population is in dataset
            if dataset.populations.contains_key(&population){
                let sequence = dataset.populations[&population].clone();
                sequences.push(sequence);
                debug!("Added population {population} to query sequences.");         
            } else {
                warn!("Population {population} was not found in the dataset.");
                continue
            }
           
        }
    }

    // ------------------------------------------------------------------------
    // input alignment was specified

    let alignment = args.input.alignment;
    if let Some( alignment ) = alignment {
        info!("Loading alignment: {:?}", alignment);
        let alignment_reader = fasta::Reader::from_file(&alignment).expect("Unable to read alignment");

        for result in alignment_reader.records() {
            let record = result.wrap_err(format!(
                "Unable to parse alignment: {:?}",
                alignment.to_str().unwrap()
            ))?;
            let sequence = Sequence::from_record(record, Some(&dataset.reference), args.mask)?;
            sequences.push(sequence);       
    
        }

    }

    debug!("Query sequences: {:?}", sequences.iter().map(|seq| seq.id.clone()).collect::<Vec<_>>());

    Ok(())
}

