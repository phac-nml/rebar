use crate::sequence::Substitution;
use crate::traits::ToYaml;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct MatchSummary {
    pub consensus_population: String,
    pub top_populations: Vec<String>,
    pub support: BTreeMap<String, usize>,
    pub private: Vec<Substitution>,
    pub conflict_ref: BTreeMap<String, usize>,
    pub conflict_alt: BTreeMap<String, usize>,
    pub total: BTreeMap<String, isize>,
}

impl ToYaml for MatchSummary {}

impl MatchSummary {
    pub fn new() -> Self {
        MatchSummary {
            consensus_population: String::new(),
            top_populations: Vec::new(),
            support: BTreeMap::new(),
            private: Vec::new(),
            conflict_ref: BTreeMap::new(),
            conflict_alt: BTreeMap::new(),
            total: BTreeMap::new(),
        }
    }
}
