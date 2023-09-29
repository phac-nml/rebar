use crate::cli;
use crate::dataset::{attributes::Name, sarscov2};
use color_eyre::eyre::{eyre, Report, Result};
use color_eyre::Help;
use log::debug;
use petgraph::dot::{Config, Dot};
use petgraph::graph::{Graph, NodeIndex};
use petgraph::visit::Dfs;
use petgraph::Direction;
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::string::ToString;

// ----------------------------------------------------------------------------
// Phylogeny Export Format

pub enum PhylogenyExportFormat {
    Dot,
    Json,
}

impl PhylogenyExportFormat {
    pub fn extension(&self) -> String {
        match self {
            PhylogenyExportFormat::Dot => String::from("dot"),
            PhylogenyExportFormat::Json => String::from("json"),
        }
    }
}

impl std::fmt::Display for PhylogenyExportFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            PhylogenyExportFormat::Dot => write!(f, "dot"),
            PhylogenyExportFormat::Json => write!(f, "json"),
        }
    }
}

pub enum PhylogenyImportFormat {
    Json,
}

// ----------------------------------------------------------------------------
// Phylogeny Import Format

impl PhylogenyImportFormat {
    pub fn extension(&self) -> String {
        match self {
            PhylogenyImportFormat::Json => String::from("json"),
        }
    }
}

impl std::fmt::Display for PhylogenyImportFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            PhylogenyImportFormat::Json => write!(f, "json"),
        }
    }
}

// ----------------------------------------------------------------------------
// Phylogeny

#[derive(Debug, Deserialize, Serialize)]
pub struct Phylogeny {
    pub graph: Graph<String, isize>,
    pub order: Vec<String>,
    pub lookup: HashMap<String, NodeIndex>,
    pub recombinants: Vec<String>,
}

impl Default for Phylogeny {
    fn default() -> Self {
        Self::new()
    }
}

impl Phylogeny {
    pub fn new() -> Self {
        Phylogeny {
            graph: Graph::new(),
            order: Vec::new(),
            lookup: HashMap::new(),
            recombinants: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.lookup.len() == 0
    }

    pub fn is_recombinant(&self, name: &str) -> Result<bool, Report> {
        let node = self.get_node(name).unwrap();
        let mut edges = self
            .graph
            .neighbors_directed(node, Direction::Incoming)
            .detach();

        // Recombinants have more than 1 incoming edge
        let mut num_edges = 0;
        while let Some(_edge) = edges.next_edge(&self.graph) {
            num_edges += 1;
        }

        if num_edges > 1 {
            Ok(true)
        } else {
            Ok(false)
        }
    }

    pub fn get_recombinants(&self) -> Result<Vec<String>, Report> {
        let mut recombinants: Vec<String> = Vec::new();

        for node in self.graph.node_indices() {
            let name = self.get_name(&node).unwrap();
            let is_recombinant = self.is_recombinant(&name).unwrap();
            if is_recombinant {
                recombinants.push(name);
            }
        }

        Ok(recombinants)
    }

    pub async fn build_graph(
        &mut self,
        args: &cli::DatasetDownloadArgs,
    ) -> Result<(), Report> {
        let (graph_data, order) = match &args.name {
            &Name::SarsCov2 => sarscov2::create_graph_data(args).await?,
            _ => {
                return Err(eyre!(
                    "Building graph for {} is not implemented.",
                    &args.name
                ))
            }
        };

        self.order = order;

        // ------------------------------------------------------------------------
        // Construct Graph

        // Add root node
        let name = "root".to_string();
        let id = self.graph.add_node(name.clone());
        self.lookup.insert(name, id);

        // todo!() Do this twice? in case lineages are accidentally out of order?

        // Add descendants
        for name in &self.order {
            let id = self.graph.add_node(name.clone());
            self.lookup.insert(name.clone(), id);
            let parents = &graph_data[&name.clone()];

            debug!("Population: {name}; Parents: {parents:?}");

            // If multiple parents add this to recombinants list
            if parents.len() > 1 {
                self.recombinants.push(name.clone());
            }
            for parent in parents {
                if !self.lookup.contains_key(parent) {
                    return Err(eyre!("Parental lineage {parent} is not in the graph.")
                        .suggestion(
                            "Are the alias_key.json and lineage_notes.txt out of sync?",
                        )
                        .suggestion("Please check if {parent} is in the alias key."));
                }
                let parent_id = self.lookup[&parent.clone()];
                self.graph.add_edge(parent_id, id, 1);
            }
        }

        Ok(())
    }

    /// import phylogeny from specified format
    pub fn import(
        dataset_dir: &Path,
        format: PhylogenyImportFormat,
    ) -> Result<Phylogeny, Report> {
        //  import path
        let mut import_path = dataset_dir.join("phylogeny");
        import_path.set_extension(format.extension());

        let phylogeny: Phylogeny = match format {
            PhylogenyImportFormat::Json => {
                let phylogeny = std::fs::read_to_string(import_path)
                    .expect("Couldn't read phylogeny {import_path:?}.");
                serde_json::from_str(&phylogeny)?
            } // ... space for more enum options for import, like Auspice perhaps
        };

        Ok(phylogeny)
    }

    /// export phylogeny to specified format
    pub fn export(
        &self,
        dataset_dir: &Path,
        format: PhylogenyExportFormat,
    ) -> Result<(), Report> {
        // output path
        let mut output_path = dataset_dir.join("phylogeny");
        output_path.set_extension(format.extension());

        debug!("Exporting phylogeny to {format}: {output_path:?}");

        // format conversion
        let output = match format {
            PhylogenyExportFormat::Dot => {
                let mut output =
                    format!("{}", Dot::with_config(&self.graph, &[Config::EdgeNoLabel]));
                // set graph id (for cytoscape)
                output = str::replace(&output, "digraph", "digraph G");
                // set horizontal (Left to Right) format for tree-like visualizer
                output =
                    str::replace(&output, "digraph {", "digraph {\n    rankdir=\"LR\";");
                output
            }
            PhylogenyExportFormat::Json => serde_json::to_string_pretty(&self)
                .unwrap_or_else(|_| panic!("Failed to export phylogeny to {format}.")),
        };

        // Write to file
        let mut file = File::create(&output_path).unwrap_or_else(|_| {
            panic!("Failed to access output phylogeny path {:?}.", &output_path)
        });
        file.write_all(output.as_bytes()).unwrap_or_else(|_| {
            panic!("Failed to write phylogeny to {:?}.", &output_path)
        });

        Ok(())
    }

    pub fn export_json(&self, dataset_dir: &Path) -> Result<(), Report> {
        let output_path = dataset_dir.join("phylogeny.json");
        let output =
            serde_json::to_string(&self).expect("Failed to export phylogeny to json.");

        let mut file = File::create(output_path).unwrap();
        file.write_all(output.as_bytes())
            .expect("Failed to export phylogeny as json.");
        Ok(())
    }

    pub fn get_descendants(&self, name: &str) -> Result<Vec<String>, Report> {
        let mut descendants = Vec::new();

        // Find the node that matches the name
        let node = self.get_node(name)?;
        // Construct a depth-first-search (Dfs)
        let mut dfs = Dfs::new(&self.graph, node);
        // Skip over self?
        // dfs.next(&self.graph);
        // Iterate over descendants
        while let Some(nx) = dfs.next(&self.graph) {
            // Get node name
            let nx_name = self.get_name(&nx).unwrap();
            descendants.push(nx_name);
        }

        Ok(descendants)
    }

    /// Get parent names of node
    pub fn get_parents(&self, name: &str) -> Result<Vec<String>, Report> {
        let mut parents = Vec::new();

        let node = self.get_node(name)?;
        let mut neighbors = self
            .graph
            .neighbors_directed(node, Direction::Incoming)
            .detach();
        while let Some(parent_node) = neighbors.next_node(&self.graph) {
            let parent_name = self.get_name(&parent_node).unwrap();
            parents.push(parent_name);
        }

        Ok(parents)
    }

    /// Get all paths from the origin node to the destination node, always traveling
    /// in the specified direction (Incoming towards root, Outgoing towards tips)/
    /// petgraph must have this already implemented, but I can't find it in docs
    pub fn get_paths(
        &self,
        origin: &str,
        dest: &str,
        direction: petgraph::Direction,
    ) -> Result<Vec<Vec<String>>, Report> {
        // container to hold the paths we've found, is a vector of vectors
        // because there might be recombinants with multiple paths
        let mut paths: Vec<Vec<String>> = Vec::new();

        // check that the origin and dest actually exist in the graph
        let origin_node = self.get_node(origin)?;
        let _dest_node = self.get_node(dest)?;

        // Check if we've reached the destination
        if origin == dest {
            paths.push(vec![origin.to_string()]);
        }
        // Otherwise, continue the search!
        else {
            let mut neighbors = self
                .graph
                .neighbors_directed(origin_node, direction)
                .detach();
            while let Some(parent_node) = neighbors.next_node(&self.graph) {
                // convert the parent graph index to a string name
                let parent_name = self.get_name(&parent_node).unwrap();

                // recursively get path of each parent to the destination
                let mut parent_paths = self.get_paths(&parent_name, dest, direction)?;

                // prepend the origin to the paths
                parent_paths
                    .iter_mut()
                    .for_each(|p| p.insert(0, origin.to_string()));

                // update the paths container to return at end of function
                for p in parent_paths {
                    paths.push(p);
                }
            }
        }

        Ok(paths)
    }

    /// NOTE: Don't think this will work with 3+ parents yet, to be tested.
    pub fn get_ancestors(&self, name: &str) -> Result<Vec<Vec<String>>, Report> {
        let mut paths = self.get_paths(name, "root", petgraph::Incoming)?;

        // remove self name (first element) from paths, and then reverse order
        // so that it's ['root'.... name]
        paths.iter_mut().for_each(|p| {
            p.remove(0);
            p.reverse();
        });

        Ok(paths)
    }

    /// Identify the most recent common ancestor shared between all node names.
    pub fn get_common_ancestor(&self, names: &Vec<String>) -> Result<String, Report> {
        // if only one node name was provided, just return it
        if names.len() == 1 {
            let common_ancestor = names[0].clone();
            return Ok(common_ancestor);
        }

        // Phase 1: Count up the ancestors shared between all named populations
        let mut ancestor_counts: HashMap<String, Vec<String>> = HashMap::new();
        let mut ancestor_depths: HashMap<String, isize> = HashMap::new();

        for name in names {
            // directly use the get_paths method over get_ancestors, because
            // get_ancestors removes the self node name from the list,
            // but some datasets have named internal nodes, so a listed
            // node could be a common ancestor!
            let ancestor_paths = self.get_paths("root", name, petgraph::Outgoing)?;

            for ancestor_path in ancestor_paths {
                for (depth, ancestor) in ancestor_path.iter().enumerate() {
                    let depth = depth as isize;
                    // add ancestor if first time encountered
                    ancestor_depths.entry(ancestor.clone()).or_insert(depth);

                    // recombinants can appear multiple times in ancestors, update
                    // depth map to use deepest one
                    if depth > ancestor_depths[ancestor] {
                        ancestor_depths.insert(ancestor.clone(), depth);
                    }
                    ancestor_counts
                        .entry(ancestor.clone())
                        .and_modify(|p| {
                            p.push(name.clone());
                            p.dedup();
                        })
                        .or_insert(vec![name.clone()]);
                }
            }
        }

        // Phase 2: Find the highest depth ancestor shared between all
        let mut common_ancestor = "root".to_string();
        let mut max_depth = 0;

        for (ancestor, populations) in ancestor_counts {
            // Which ancestors were found in all populations?
            if populations.len() == names.len() {
                // Which ancestor has the max depth?
                let depth = ancestor_depths[&ancestor];
                if depth > max_depth {
                    max_depth = depth;
                    common_ancestor = ancestor;
                }
            }
        }

        Ok(common_ancestor)
    }

    /// Identify the most recent ancestor that is a recombinant.
    pub fn get_recombinant_ancestor(&self, name: &str) -> Result<Option<String>, Report> {
        let mut recombinant: Option<String> = None;

        let ancestor_paths = self.get_paths(name, "root", petgraph::Incoming)?;

        for path in ancestor_paths {
            for name in path {
                if self.recombinants.contains(&name) {
                    recombinant = Some(name.to_string());
                    break;
                }
            }
            if recombinant.is_some() {
                break;
            }
        }

        Ok(recombinant)
    }

    pub fn get_node(&self, name: &str) -> Result<NodeIndex, Report> {
        if self.lookup.contains_key(name) {
            let node = self.lookup[name];
            Ok(node)
        } else {
            Err(eyre!("Node {name} is not found in the phylogeny."))
        }
    }

    pub fn get_name(&self, node: &NodeIndex) -> Option<String> {
        for (name, node_l) in &self.lookup {
            if node == node_l {
                return Some(name.clone());
            }
        }

        None
    }
}
