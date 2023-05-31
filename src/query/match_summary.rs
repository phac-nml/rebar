use crate::sequence::Substitution;
use crate::traits::ToYaml;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct MatchSummary {
    pub consensus_population: String,
    pub top_populations: Vec<String>,
    pub substitutions: Vec<Substitution>,
    pub support: BTreeMap<String, Vec<Substitution>>,
    pub private: Vec<Substitution>,
    pub conflict_ref: BTreeMap<String, Vec<Substitution>>,
    pub conflict_alt: BTreeMap<String, Vec<Substitution>>,
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
        for (pop, _count) in &total_order {
            let subs = &self.support[*pop];
            let count = subs.len();
            support_order.push(format!(
                "{}:\n    - count: {}\n    - substitutions: {}",
                pop,
                count,
                subs.iter().join(", ")
            ));
        }

        // conflict_ref
        let mut conflict_ref_order: Vec<String> = Vec::new();
        for (pop, _count) in &total_order {
            let subs = &self.conflict_ref[*pop];
            let count = subs.len();
            conflict_ref_order.push(format!(
                "{}:\n    - count: {}\n    - substitutions: {}",
                pop,
                count,
                subs.iter().join(", ")
            ));
        }

        // conflict_alt
        let mut conflict_alt_order: Vec<String> = Vec::new();
        for (pop, _count) in &total_order {
            let subs = &self.conflict_alt[*pop];
            let count = subs.len();
            conflict_alt_order.push(format!(
                "{}:\n    - count: {}\n    - substitutions: {}",
                pop,
                count,
                subs.iter().join(", ")
            ));
        }

        // Pretty string formatting for yaml
        let total_order = total_order
            .iter()
            .map(|(pop, count)| format!("{}:\n    - count: {}", &pop, &count))
            .collect::<Vec<_>>();

        // private count and list
        let private_order = format!(
            "  - count: {}\n    - substitutions: {}",
            self.private.len(),
            self.private.iter().join(", ")
        );

        format!(
            "consensus_population: {}
top_populations:\n  - {}
recombinant: {}
substitutions: {}
total:\n  {}
support:\n  {}
conflict_ref:\n  {}
conflict_alt:\n  {}
private:\n  {}",
            self.consensus_population,
            self.top_populations.join("\n  - "),
            self.recombinant.clone().unwrap_or("None".to_string()),
            self.substitutions.iter().join(", "),
            total_order.join("\n  "),
            support_order.join("\n  "),
            conflict_ref_order.join("\n  "),
            conflict_alt_order.join("\n  "),
            private_order,
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
            substitutions: Vec::new(),
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
