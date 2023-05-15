use petgraph::graph::{Graph, NodeIndex};
use std::collections::BTreeMap;
use std::default::Default;

#[derive(Debug, Default)]
pub struct Phylogeny {
    pub graph: Graph<String, isize>,
    pub lookup: BTreeMap<String, NodeIndex>,
}
