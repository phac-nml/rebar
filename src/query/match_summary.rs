use crate::sequence::Substitution;
use crate::traits::ToYaml;
use itertools::Itertools;
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
    pub recombinant: Option<String>,
}

impl ToYaml for MatchSummary {
    fn to_yaml(&self) -> String {
        // Order the population lists from 'best' to 'worst'

        // total
        let mut total_order = self.total.iter().collect::<Vec<(_, _)>>();
        total_order.sort_by(|a, b| b.1.cmp(a.1));

        // support
        let mut support_order: Vec<String> = Vec::new();
        for (pop, _total) in &total_order {
            let count = self.support[*pop];
            support_order.push(format!("{}: {}", pop, count))
        }

        // conflict_ref
        let mut conflict_ref_order: Vec<String> = Vec::new();
        for (pop, _total) in &total_order {
            let count = self.conflict_ref[*pop];
            conflict_ref_order.push(format!("{}: {}", pop, count))
        }

        // conflict_alt
        let mut conflict_alt_order: Vec<String> = Vec::new();
        for (pop, _total) in &total_order {
            let count = self.conflict_alt[*pop];
            conflict_alt_order.push(format!("{}: {}", pop, count))
        }

        // Pretty string formatting for yaml
        let total_order = total_order
            .iter()
            .map(|(pop, count)| format!("{}: {}", &pop, &count))
            .collect::<Vec<_>>();

        format!(
            "consensus_population: {}
top_populations: \n- {}
recombinant: {}
total:\n  {}
support:\n  {}
conflict_ref:\n  {}
conflict_alt:\n  {}
private:\n  {}",
            self.consensus_population,
            self.top_populations.join("\n- "),
            self.recombinant.clone().unwrap_or("None".to_string()),
            total_order.join("\n  "),
            support_order.join("\n  "),
            conflict_ref_order.join("\n  "),
            conflict_alt_order.join("\n  "),
            self.private.iter().join("\n- ")
        )
    }
}

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
            recombinant: None,
        }
    }
}

impl Default for MatchSummary {
    fn default() -> Self {
        Self::new()
    }
}
