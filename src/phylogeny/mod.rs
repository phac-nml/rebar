use crate::utils;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use color_eyre::Help;
use itertools::Itertools;
use log::debug;
use petgraph::dot::{Config, Dot};
use petgraph::graph::{Graph, NodeIndex};
use petgraph::visit::{Dfs, IntoNodeReferences};
use petgraph::Direction;
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::string::ToString;

// ----------------------------------------------------------------------------
// Phylogeny

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Phylogeny {
    pub graph: Graph<String, isize>,
    pub order: Vec<String>,
    pub recombinants: Vec<String>,
    pub recombinants_all: Vec<String>,
    pub non_recombinants_all: Vec<String>,
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
            recombinants: Vec::new(),
            recombinants_all: Vec::new(),
            non_recombinants_all: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.graph.node_count() == 0
    }

    /// Return true if a node name is a recombinant.
    pub fn is_recombinant(&self, name: &str) -> Result<bool, Report> {
        let node = self.get_node(name)?;
        let mut edges = self.graph.neighbors_directed(node, Direction::Incoming).detach();

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

    /// Get recombinant node names.
    pub fn get_recombinants(&self) -> Result<Vec<String>, Report> {
        let mut recombinants: Vec<String> = Vec::new();

        for node in self.graph.node_indices() {
            let name = self.get_name(&node)?;
            let is_recombinant = self.is_recombinant(&name)?;
            if is_recombinant {
                recombinants.push(name);
            }
        }

        Ok(recombinants)
    }

    /// Get recombinant node names and all descendants.
    pub fn get_recombinants_all(&self) -> Result<Vec<String>, Report> {
        let mut recombinants: Vec<String> = Vec::new();

        for node in self.graph.node_indices() {
            let name = self.get_name(&node)?;
            let result = self.get_recombinant_ancestor(&name)?;
            if result.is_some() {
                recombinants.push(name)
            }
        }

        Ok(recombinants)
    }

    /// Get non-recombinants
    pub fn get_non_recombinants_all(&self) -> Result<Vec<String>, Report> {
        let mut non_recombinants: Vec<String> = Vec::new();

        for node in self.graph.node_indices() {
            let name = self.get_name(&node)?;
            let result = self.get_recombinant_ancestor(&name)?;
            if result.is_none() {
                non_recombinants.push(name)
            }
        }

        Ok(non_recombinants)
    }

    /// Prune a clade from the graph.
    pub fn prune(&self, name: &str) -> Result<Phylogeny, Report> {
        let mut phylogeny = self.clone();

        // Find the node that matches the name
        let descendants = self.get_descendants(name)?;

        for name in &descendants {
            debug!("Pruning node: {name}");
            let node = self.get_node(name)?;
            // Remove from graph
            let removal = phylogeny.graph.remove_node(node);
            // If we have already checked that the node exists, why does this fail?
            // Ex. XBB knockout, fails on XBL
            if removal.is_none() {
                debug!("Node {node:?} was not found in the graph.");
            }
        }

        // Remove from order
        phylogeny.order.retain(|n| !descendants.contains(n));
        // Remove from recombinants
        phylogeny.recombinants.retain(|r| !descendants.contains(r));

        Ok(phylogeny)
    }

    /// Read phylogeny from file.
    pub fn read(path: &Path) -> Result<Phylogeny, Report> {
        let phylogeny = std::fs::read_to_string(path)
            .wrap_err_with(|| "Failed to read file: {path:?}.")?;
        let phylogeny = serde_json::from_str(&phylogeny)
            .wrap_err_with(|| "Failed to parse file: {path:?}.")?;

        Ok(phylogeny)
    }

    /// Write phylogeny to file.
    pub fn write(&self, output_path: &Path) -> Result<(), Report> {
        // Create output file
        let mut file = File::create(output_path)?;

        // .unwrap_or_else(|_| {
        //     return Err(eyre!("Failed to create file: {:?}.", &output_path))
        // });
        let ext = utils::path_to_ext(Path::new(output_path))?;

        // format conversion
        let output = match ext.as_str() {
            // ----------------------------------------------------------------
            // DOT file for graphviz
            "dot" => {
                let mut output =
                    format!("{}", Dot::with_config(&self.graph, &[Config::EdgeNoLabel]));
                // set graph id (for cytoscape)
                output = str::replace(&output, "digraph", "digraph G");
                // set horizontal (Left to Right) format for tree-like visualizer
                output =
                    str::replace(&output, "digraph {", "digraph {\n    rankdir=\"LR\";");
                output
            }
            // ----------------------------------------------------------------
            // JSON for rebar
            "json" => serde_json::to_string_pretty(&self)
                .unwrap_or_else(|_| panic!("Failed to parse: {self:?}")),
            _ => {
                return Err(eyre!(
                    "Phylogeny write for extension .{ext} is not supported."
                )
                .suggestion("Please try .json or .dot instead."))
            }
        };

        // Write to file
        file.write_all(output.as_bytes())
            .unwrap_or_else(|_| panic!("Failed to write file: {:?}.", &output_path));

        Ok(())
    }

    // Reminder, this function will also include name (the parent)
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
            let nx_name = self.get_name(&nx)?;
            descendants.push(nx_name);
        }

        Ok(descendants)
    }

    /// Get parent names of node
    pub fn get_parents(&self, name: &str) -> Result<Vec<String>, Report> {
        let mut parents = Vec::new();

        let node = self.get_node(name)?;
        let mut neighbors =
            self.graph.neighbors_directed(node, Direction::Incoming).detach();
        while let Some(parent_node) = neighbors.next_node(&self.graph) {
            let parent_name = self.get_name(&parent_node)?;
            parents.push(parent_name);
        }

        Ok(parents)
    }

    /// Get problematic recombinants, where the parents are not sister taxa.
    /// They might be parent-child instead.
    pub fn get_problematic_recombinants(&self) -> Result<Vec<String>, Report> {
        let mut problematic_recombinants = Vec::new();

        for recombinant in &self.recombinants {
            let parents = self.get_parents(recombinant)?;
            for i1 in 0..parents.len() - 1 {
                let p1 = &parents[i1];
                for p2 in parents.iter().skip(i1 + 1) {
                    let mut descendants = self.get_descendants(p2)?;
                    let ancestors =
                        self.get_ancestors(p2)?.into_iter().flatten().collect_vec();
                    descendants.extend(ancestors);

                    if descendants.contains(p1) {
                        problematic_recombinants.push(recombinant.clone());
                        break;
                    }
                }
            }
        }

        Ok(problematic_recombinants)
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
            let mut neighbors =
                self.graph.neighbors_directed(origin_node, direction).detach();
            while let Some(parent_node) = neighbors.next_node(&self.graph) {
                // convert the parent graph index to a string name
                let parent_name = self.get_name(&parent_node)?;

                // recursively get path of each parent to the destination
                let mut parent_paths = self.get_paths(&parent_name, dest, direction)?;

                // prepend the origin to the paths
                parent_paths.iter_mut().for_each(|p| p.insert(0, origin.to_string()));

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
        for (idx, n) in self.graph.node_references() {
            if n == name {
                return Ok(idx);
            }
        }
        Err(eyre!("Name {name} is not in the phylogeny."))
    }

    pub fn get_name(&self, node: &NodeIndex) -> Result<String, Report> {
        for (idx, n) in self.graph.node_references() {
            if &idx == node {
                return Ok(n.clone());
            }
        }
        Err(eyre!("Node {node:?} is not in the phylogeny."))
    }
}
