use crate::traits::ToYaml;
use serde::{Deserialize, Serialize};
use std::default::Default;

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Sequence {
    pub id: String,
    pub seq: Vec<char>,
}

impl ToYaml for Sequence {}

impl Sequence {
    pub fn new() -> Self {
        Sequence {
            id: String::new(),
            seq: Vec::new(),
        }
    }

    pub fn from_record(record: bio::io::fasta::Record) -> Self {
        Sequence {
            id: record.id().to_string(),
            seq: record.seq().iter().map(|b| *b as char).collect(),
        }
    }
}
