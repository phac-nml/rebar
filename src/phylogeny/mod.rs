use crate::dataset;
use crate::utils;
use color_eyre::Report;
use csv;
use itertools::Itertools;
use log::{debug, info};
use petgraph::dot::{Config, Dot};
use petgraph::graph::{Graph, NodeIndex};
use petgraph::visit::Dfs;
use petgraph::Direction;
//use serde;
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::string::ToString;

#[derive(serde::Deserialize, Debug)]
struct LineageNotesRow<'a> {
    lineage: &'a str,
    _description: &'a str,
}

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

    pub fn is_recombinant(&self, name: &String) -> Result<bool, Report> {
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
        dataset_name: &dataset::Name,
        dataset_dir: &Path,
    ) -> Result<(), Report> {
        let (graph_data, order) = match dataset_name {
            &dataset::Name::SarsCov2 => create_sarscov2_graph_data(dataset_dir).await?,
            _ => return Ok(()),
        };

        self.order = order;

        // ------------------------------------------------------------------------
        // Construct Graph

        // Add root node
        let name = "root".to_string();
        let id = self.graph.add_node(name.clone());
        self.lookup.insert(name, id);

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

        info!("Exporting phylogeny to {format}: {output_path:?}");

        // format conversion
        let output = match format {
            PhylogenyExportFormat::Dot => {
                let mut output =
                    format!("{}", Dot::with_config(&self.graph, &[Config::EdgeNoLabel]));
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

    pub fn get_descendants(&self, name: &String) -> Result<Vec<String>, Report> {
        let mut descendants = Vec::new();

        // Find the node that matches the name
        let node = self
            .get_node(name)
            .unwrap_or_else(|| panic!("Couldn't find node name in phylogeny: {}", name));
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

    pub fn get_ancestors(&self, name: &String) -> Result<Vec<Vec<String>>, Report> {
        let mut ancestors = Vec::new();

        // Copy graph so we can mutate it here with reverse
        let mut graph = self.graph.clone();

        // Construct a backwards depth-first-search (Dfs)
        graph.reverse();
        let node = self.get_node(name).unwrap();
        let mut dfs = Dfs::new(&graph, node);

        // Walk to the root, there might be multiple paths (recombinants)
        let mut path_nodes = Vec::new();
        let mut prev_name = String::new();

        // Skip self?
        //dfs.next(&graph);
        while let Some(nx) = dfs.next(&graph) {
            // Get node name
            let nx_name = self.get_name(&nx).unwrap();

            // If the previous node name was root, that means we topped
            // out the search in the last iter, but still have alternate
            // recombinant paths to deal with
            if prev_name == "root" {
                // Add the topped out path to our list of paths
                path_nodes.reverse();
                ancestors.push(path_nodes.clone());

                // Initialize vector for new paths
                //path_nodes = Vec::new();

                // Recursive search, swap graph back and forth
                graph.reverse();
                let nx_ancestors = self.get_ancestors(&nx_name).unwrap();
                for ancestor_nodes in &nx_ancestors {
                    ancestors.push(ancestor_nodes.clone());
                }
                graph.reverse();
            }

            path_nodes.push(nx_name.clone());

            prev_name = nx_name;
        }

        if prev_name == "root" {
            path_nodes.reverse();
            ancestors.push(path_nodes.clone());
        }

        // Restore original graph order
        graph.reverse();
        Ok(ancestors)
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
            let ancestor_paths = self.get_ancestors(name).unwrap();
            for ancestor_path in ancestor_paths {
                for (depth, ancestor) in ancestor_path.iter().enumerate() {
                    let depth = depth as isize;
                    ancestor_depths.entry(ancestor.clone()).or_insert(depth);
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

    pub fn get_node(&self, name: &String) -> Option<NodeIndex> {
        if self.lookup.contains_key(name) {
            let node = self.lookup[name];
            return Some(node);
        }

        None
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

pub async fn download_lineage_notes(dataset_dir: &Path) -> Result<PathBuf, Report> {
    let decompress = false;
    let url = dataset::SARSCOV2_LINEAGE_NOTES_URL;
    let output_path = dataset_dir.join("lineage_notes.txt");
    utils::download_file(url, &output_path, decompress).await?;

    Ok(output_path)
}

pub async fn download_alias_key(dataset_dir: &Path) -> Result<PathBuf, Report> {
    // https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json

    let decompress = false;
    let url = dataset::SARSCOV2_ALIAS_KEY_URL;
    let output_path = dataset_dir.join("alias_key.json");
    utils::download_file(url, &output_path, decompress).await?;

    Ok(output_path)
}

pub async fn import_alias_key(
    dataset_dir: &Path,
) -> Result<HashMap<String, Vec<String>>, Report> {
    // download
    let alias_key_path = download_alias_key(dataset_dir)
        .await
        .expect("Couldn't download alias_key from url.");

    // read to json
    let alias_key_str =
        std::fs::read_to_string(alias_key_path).expect("Couldn't read alias_key file.");
    let alias_key_val: serde_json::Value =
        serde_json::from_str(&alias_key_str).expect("Couldn't convert alias_key to json");
    // deserialize to object (raw, mixed types)
    let alias_key_raw: serde_json::Map<String, serde_json::Value> = alias_key_val
        .as_object()
        .expect("Couldn't convert alias_key json to json Map")
        .clone();

    // This should probably be a custom deserializer fn for brevity,
    //   but I don't know how to do that yet :)
    let mut alias_key: HashMap<String, Vec<String>> = HashMap::new();

    for (alias, lineage) in &alias_key_raw {
        let mut lineage_paths: Vec<String> = Vec::new();

        // Consistify the alias key types
        match lineage.as_array() {
            // If array, this is a recombinant alias with multiple parents.
            Some(parents) => {
                for parent in parents {
                    let parent =
                        parent.as_str().expect("Couldn't convert parent to str.");
                    // Strip the wildcard asterisks from lineage name
                    let parent_clean = str::replace(parent, "*", "");
                    lineage_paths.push(parent_clean);
                }
            }
            // Otherwise, it might be a string
            None => {
                let mut lineage_path = lineage
                    .as_str()
                    .expect("Couldn't convert lineage to str.")
                    .to_string();
                // If there is not lineage_path (ex. "" for A, B), set to self
                if lineage_path.is_empty() {
                    lineage_path = alias.clone();
                }
                lineage_paths.push(lineage_path);
            }
        }

        alias_key.insert(alias.clone(), lineage_paths);
    }

    Ok(alias_key)
}

pub async fn create_sarscov2_graph_data(
    dataset_dir: &Path,
) -> Result<(HashMap<String, Vec<String>>, Vec<String>), Report> {
    // ------------------------------------------------------------------------
    // Download and import data

    // Import the alias key
    let alias_key = import_alias_key(dataset_dir).await?;

    // Import the lineage notes
    let lineage_notes_path = download_lineage_notes(dataset_dir)
        .await
        .expect("Couldn't download lineage notes from url.");
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(lineage_notes_path)?;
    let _headers = reader.headers()?;

    // ------------------------------------------------------------------------
    // Map parent-child relationships

    let mut graph_data: HashMap<String, Vec<String>> = HashMap::new();
    let mut graph_order: Vec<String> = Vec::new();

    for result in reader.records() {
        let record = result?;

        let row: LineageNotesRow = record.deserialize(None)?;
        let lineage = row.lineage.to_string();

        // Lineages that start with '*' have been withdrawn
        if lineage.starts_with('*') {
            continue;
        }

        let parents = get_lineage_parents(lineage.clone(), &alias_key).unwrap();

        //debug!("Lineage: {lineage}; Parents: {parents:?}");

        graph_order.push(lineage.to_string().clone());
        graph_data.insert(lineage.clone(), parents.clone());
    }

    Ok((graph_data, graph_order))
}

pub fn get_lineage_parents(
    lineage: String,
    alias_key: &HashMap<String, Vec<String>>,
) -> Result<Vec<String>, Report> {
    let mut parents: Vec<String> = Vec::new();

    // If Recombinant with multiple parents, if so, it will be in the alias
    // key with the parents listed.
    if alias_key.contains_key(&lineage) {
        let lineage_paths = &alias_key[&lineage];
        if lineage_paths.len() > 1 {
            // Dedup in case multiple breakpoints/parents
            parents = lineage_paths.clone().into_iter().unique().collect();
            return Ok(parents);
        }
    }

    // Otherwise, single parent
    let decompress = decompress_lineage(&lineage, alias_key).unwrap();

    // Ex. BA.5.2 -> ["BA", "5", "2"]
    let decompress_parts = decompress
        .split('.')
        .map(|p| p.to_string())
        .collect::<Vec<_>>();

    // If just 1 part, parent is root (ex. A)
    let mut parent = String::from("root");
    if decompress_parts.len() > 1 {
        parent = decompress_parts[0..(decompress_parts.len() - 1)].join(".");
    }

    // Compress the full parent back down with aliases
    parent = compress_lineage(&parent, alias_key).unwrap();
    parents.push(parent);

    Ok(parents)
}

pub fn compress_lineage(
    lineage: &String,
    alias_key: &HashMap<String, Vec<String>>,
) -> Result<String, Report> {
    // By default, set compression level to self
    let mut compress = lineage.to_string();

    // Reverse the alias-> lineage path lookup
    let mut alias_key_rev: HashMap<String, String> = HashMap::new();

    for (alias, lineage_paths) in alias_key {
        // Skip over recombinants with multiple parents, don't need their lookup
        if lineage_paths.len() > 1 {
            continue;
        }
        let lineage_path = lineage_paths[0].clone();
        alias_key_rev.insert(lineage_path, alias.clone());
    }

    // Ex. BA.5.2 -> ["BA", "5", "2"]
    let compress_parts = compress
        .split('.')
        .map(|p| p.to_string())
        .collect::<Vec<_>>();

    if compress_parts.len() > 1 {
        for i in (0..compress_parts.len()).rev() {
            let compress_subset = compress_parts[0..i].join(".");

            if alias_key_rev.contains_key(&compress_subset) {
                compress = alias_key_rev[&compress_subset].clone();
                // Get the suffix that was chopped off in subset
                let compress_suffix = &compress_parts[i..];

                // Add the suffix
                if !compress_suffix.is_empty() {
                    compress = format!["{compress}.{}", compress_suffix.join(".")];
                }
                break;
            }
        }
    }

    Ok(compress)
}

/// Decompresses a SARS-CoV-2 lineage name.
///
/// Convert a compressed, SARS-CoV-2 lineage name into it's full unaliased form. q
///
/// # Arguments
///
///  * `lineage` | `String` | A String of a SARS-CoV-2 lineage name.
///  * `alias_key` | `&HashMap<String, Vec<String>>` | A mapping of SARS-CoV-2 aliases to decompressed lineage paths.
///
/// # Basic usage:
///
/// ```
/// decompress_lineage("BA.5.2", alias_key).unwrap();
/// ```
///
/// # Example
///
/// ```no run
/// let alias_key = import_alias_key(".").unwrap();
/// decompress_lineage("BA.5.2", alias_key).unwrap();
/// ```
pub fn decompress_lineage(
    lineage: &str,
    alias_key: &HashMap<String, Vec<String>>,
) -> Result<String, Report> {
    // By default, set full path to lineage
    let mut decompress = lineage.to_string();

    // Split lineage into levels, Ex. BA.5.2 = ["BA", "5", "2"]
    // Can be a maximum of 4 levels before aliasing
    let mut lineage_level = 0;
    let mut lineage_parts = vec![String::new(); 4];

    for (i, level) in lineage.split('.').enumerate() {
        lineage_parts[i] = level.to_string();
        lineage_level = i + 1;
    }

    // Get the first letter prefix, Ex. BA.5.2 = "BA"
    let lineage_prefix = lineage_parts[0].clone();
    // If there were multiple parts, get suffix (Ex. "BA" and "5.2")
    let lineage_suffix = lineage_parts[1..lineage_level].join(".");

    // Decompressing logic
    if alias_key.contains_key(&lineage_prefix) {
        let lineage_paths = &alias_key[&lineage_prefix];
        // Not multiple recombinant parents
        if lineage_paths.len() == 1 {
            decompress = lineage_paths[0].clone();
            // Add back our suffix numbers
            if lineage_level > 1 {
                decompress = format!("{decompress}.{lineage_suffix}");
            }
        }
    }

    Ok(decompress)
}
