use crate::phylogeny::Phylogeny;
use color_eyre::eyre::{Report, Result};
use itertools::Itertools;
use std::path::PathBuf;

pub fn build() -> Result<Phylogeny, Report> {
    let mut phylogeny = Phylogeny::new();
    phylogeny.order = vec!["A", "B", "C"].into_iter().map(String::from).collect_vec();

    // Add root node
    let name = "root".to_string();
    phylogeny.graph.add_node(name.clone());
    let root_id = phylogeny.get_node(&name)?;

    let name = "A".to_string();
    phylogeny.graph.add_node(name.clone());
    let a_id = phylogeny.get_node(&name)?;
    phylogeny.graph.add_edge(root_id, a_id, 1);

    let name = "B".to_string();
    phylogeny.graph.add_node(name.clone());
    let b_id = phylogeny.get_node(&name)?;
    phylogeny.graph.add_edge(root_id, b_id, 1);

    let name = "C".to_string();
    phylogeny.graph.add_node(name.clone());
    let c_id = phylogeny.get_node(&name)?;
    phylogeny.graph.add_edge(a_id, c_id, 1);
    phylogeny.graph.add_edge(b_id, c_id, 1);
    phylogeny.recombinants.push(name.clone());

    phylogeny.recombinants_all = phylogeny.get_recombinants_all()?;

    let output_path = PathBuf::from("phylogeny.json");
    phylogeny.write(&output_path)?;

    Ok(phylogeny)
}
