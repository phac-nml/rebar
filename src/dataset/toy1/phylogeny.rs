use crate::phylogeny::Phylogeny;
use color_eyre::eyre::{Report, Result};

pub fn build() -> Result<Phylogeny, Report> {
    let mut phylogeny = Phylogeny::new();

    // Add root node
    let name = "root".to_string();
    let root_id = phylogeny.graph.add_node(name.clone());

    // Add A node
    let name = "A".to_string();
    let a_id = phylogeny.graph.add_node(name.clone());
    phylogeny.graph.add_edge(root_id, a_id, 1);

    // Add B node
    let name = "B".to_string();
    let b_id = phylogeny.graph.add_node(name.clone());
    phylogeny.graph.add_edge(root_id, b_id, 1);

    // Add C node
    let name = "C".to_string();
    let c_id = phylogeny.graph.add_node(name.clone());
    phylogeny.graph.add_edge(root_id, c_id, 1);

    // Add recombinant D node
    let name = "D".to_string();
    let d_id = phylogeny.graph.add_node(name.clone());
    phylogeny.graph.add_edge(a_id, d_id, 1);
    phylogeny.graph.add_edge(b_id, d_id, 1);

    // Add recursive recombinant E node
    let name = "E".to_string();
    let e_id = phylogeny.graph.add_node(name.clone());
    phylogeny.graph.add_edge(d_id, e_id, 1);
    phylogeny.graph.add_edge(c_id, e_id, 1);

    Ok(phylogeny)
}
