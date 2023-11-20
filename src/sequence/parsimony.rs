use crate::sequence::{Sequence, Substitution};
use color_eyre::eyre::{Report, Result};
use indoc::formatdoc;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

// ----------------------------------------------------------------------------
// Population Parsimony Summary

/// Summarize support and conflicts between two sequences.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Summary {
    pub support: Vec<Substitution>,
    pub conflict_ref: Vec<Substitution>,
    pub conflict_alt: Vec<Substitution>,
    pub score: isize,
}

impl Summary {
    pub fn new() -> Self {
        Summary {
            support: Vec::new(),
            conflict_ref: Vec::new(),
            conflict_alt: Vec::new(),
            score: 0,
        }
    }

    /// Summarize support and conflicts between two sequences.
    pub fn from_sequence(
        sequence: &Sequence,
        query: &Sequence,
        coordinates: Option<&Vec<usize>>,
    ) -> Result<Self, Report> {
        let mut parsimony_summary = Summary::new();

        // get all the substitutions found in this query
        // exclude missing and deletion coordinates
        let mut query_subs = query
            .substitutions
            .iter()
            .filter(|sub| {
                !sequence.missing.contains(&sub.coord)
                    && !sequence.deletions.contains(&sub.to_deletion())
            })
            .collect_vec();

        // get all the substitutions found in the sequence
        let mut seq_subs = sequence.substitutions.clone();

        // optionally filter coordinates
        if let Some(coordinates) = coordinates {
            query_subs.retain(|sub| coordinates.contains(&sub.coord));
            seq_subs.retain(|sub| coordinates.contains(&sub.coord));
        }

        // support: sub in seq that is also in query
        // conflict_alt: sub in eq that is not in query
        seq_subs.iter().for_each(|sub| {
            if query_subs.contains(&sub) {
                parsimony_summary.support.push(*sub);
            } else {
                parsimony_summary.conflict_alt.push(*sub);
            }
        });

        // conflict_ref: sub in query that is not in seq
        parsimony_summary.conflict_ref = query_subs
            .into_iter()
            .filter(|sub| !seq_subs.contains(sub))
            .cloned()
            .collect_vec();

        // score: support - conflict_alt - conflict_ref
        // why did we previously use only conflict_ref and not conflict_alt?
        // individual isize conversion otherwise: "attempt to subtract with overflow"
        parsimony_summary.score = parsimony_summary.support.len() as isize
            - parsimony_summary.conflict_ref.len() as isize
            - parsimony_summary.conflict_alt.len() as isize;

        Ok(parsimony_summary)
    }

    pub fn pretty_print(&self) -> String {
        formatdoc!(
            "score:\n  {}
            support:\n  {}
            conflict_ref:\n  {}
            conflict_alt:\n  {}",
            self.score,
            self.support.iter().join(", "),
            self.conflict_ref.iter().join(", "),
            self.conflict_alt.iter().join(", "),
        )
    }
}

impl Default for Summary {
    fn default() -> Self {
        Self::new()
    }
}
