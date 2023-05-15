use petgraph::graph::{Graph, NodeIndex};
use std::default::Default;
use std::collections::BTreeMap;

#[derive(serde::Deserialize, Debug)]
struct LineageNotesRow<'a> {
    lineage: &'a str,
    _description: &'a str,
}

#[derive(Debug, Default)]
pub struct Phylogeny {
    pub graph: Graph<String, isize>,
    pub lookup: BTreeMap<String, NodeIndex>,
}