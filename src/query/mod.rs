pub mod match_summary;

use crate::dataset::Dataset;
use crate::sequence::Sequence;
use bio::io::fasta;
use color_eyre::eyre::{Report, Result, WrapErr};
use log::info;
use std::collections::BTreeMap;
use std::path::PathBuf;

#[allow(dead_code)]
pub struct Query {
    pub sequences: BTreeMap<String, Sequence>,
}

impl Query {
    pub fn load(
        alignment: PathBuf,
        dataset: &Dataset,
        mask: usize,
    ) -> Result<Query, Report> {
        info!("Loading query alignment: {:?}", alignment);
        let alignment_reader =
            fasta::Reader::from_file(&alignment).expect("Unable to read alignment");

        let mut sequences = BTreeMap::new();

        for result in alignment_reader.records() {
            let record = result.wrap_err(format!(
                "Unable to parse alignment: {:?}",
                alignment.to_str().unwrap()
            ))?;
            let sequence = Sequence::from_record(record, Some(&dataset.reference), mask)?;
            sequences.insert(sequence.id.clone(), sequence);
        }

        let query = Query { sequences };

        Ok(query)
    }
}
