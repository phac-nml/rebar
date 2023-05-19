use crate::dataset::Dataset;
use crate::sequence::Sequence;
use bio::io::fasta;
use color_eyre::eyre::Report;
use log::info;
use std::collections::BTreeMap;
use std::path::PathBuf;

#[allow(dead_code)]
pub struct Query {
    sequences: BTreeMap<String, Sequence>,
}

impl Query {
    pub fn load(
        alignment: PathBuf,
        dataset: &Dataset,
        mask: usize,
    ) -> Result<Query, Report> {
        info!("Loading query alignment: {:?}", alignment);
        let alignment =
            fasta::Reader::from_file(alignment).expect("Unable to load alignment");

        let mut sequences = BTreeMap::new();

        for result in alignment.records() {
            let record = result?;
            let sequence = Sequence::from_record(record, Some(&dataset.reference), mask)?;
            sequences.insert(sequence.id.clone(), sequence);
        }

        let query = Query { sequences };

        Ok(query)
    }
}
