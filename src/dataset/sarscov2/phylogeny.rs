use crate::{dataset, phylogeny::Phylogeny, utils::table::Table};
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use color_eyre::Help;
use itertools::Itertools;
use log::{debug, info, warn};
use std::collections::BTreeMap;
use std::path::Path;

pub async fn build(
    summary: &mut dataset::attributes::Summary,
    output_dir: &Path,
) -> Result<Phylogeny, Report> {
    // ------------------------------------------------------------------------
    // Download
    // ------------------------------------------------------------------------

    // ------------------------------------------------------------------------
    // Lineage Notes
    let output_path = output_dir.join("lineage_notes.txt");
    info!("Downloading lineage notes: {output_path:?}");

    let remote_file = if summary.misc.contains_key("lineage_notes") {
        dataset::download::snapshot(&summary.misc["lineage_notes"], &output_path).await?
    } else {
        dataset::sarscov2::download::lineage_notes(&summary.tag, &output_path).await?
    };
    summary.misc.insert("lineage_notes".to_string(), remote_file);

    // ------------------------------------------------------------------------
    // Alias key
    let output_path = output_dir.join("alias_key.json");
    info!("Downloading alias key: {output_path:?}");

    let remote_file = if summary.misc.contains_key("alias_key") {
        dataset::download::snapshot(&summary.misc["alias_key"], &output_path).await?
    } else {
        dataset::sarscov2::download::alias_key(&summary.tag, &output_path).await?
    };
    summary.misc.insert("alias_key".to_string(), remote_file);

    // ------------------------------------------------------------------------
    // Alias
    // ------------------------------------------------------------------------

    // read alias key into Map
    let alias_key_path = &summary.misc["alias_key"].local_path;
    let alias_key_file_name = alias_key_path.file_name().unwrap().to_str().unwrap();
    let mut alias_key = read_alias_key(alias_key_path)?;

    // MANUAL!
    // These are potential conflicts in the prior information, that will
    // cause validation to fail

    // Potentially, XAS needs deletions to resolve between BA.4/BA.5 parent?
    warn!("Changing XAS designated parent from BA.5 to BA.4: https://github.com/cov-lineages/pango-designation/issues/882");
    let recombinant = "XAS".to_string();
    let parents = vec!["BA.4".to_string(), "BA.2".to_string()];
    alias_key.insert(recombinant, parents);

    warn!("Relaxing XBB designated parent from BM.1.1.1 to BA.2.75.3");
    let recombinant = "XBB".to_string();
    let parents = vec!["BJ.1".to_string(), "BA.2.75.3".to_string()];
    alias_key.insert(recombinant, parents);

    warn!("Relaxing XBF designated parent from BA.5.2.3 to BA.5.2");
    let recombinant = "XBF".to_string();
    let parents = vec!["BA.5.2".to_string(), "CJ.1".to_string()];
    alias_key.insert(recombinant, parents);

    // ------------------------------------------------------------------------
    // Lineage Notes
    // ------------------------------------------------------------------------

    // read lineage notes into Table
    let lineage_notes_path = &summary.misc["lineage_notes"].local_path;
    let lineage_notes_file_name =
        lineage_notes_path.file_name().unwrap().to_str().unwrap();
    let lineage_table = Table::read(lineage_notes_path)?;
    // identify which column is 'Lineage'
    let lineage_col_i = lineage_table.header_position("Lineage")?;
    // get a list of lineages in the notes
    let notes_lineages = lineage_table
        .rows
        .iter()
        .filter(|row| {
            !row[lineage_col_i].starts_with('*') && row[lineage_col_i].is_empty()
        })
        .map(|row| row[lineage_col_i].to_string())
        .collect_vec();

    // read populations fasta, to check if any lineages are missing in notes
    let populations_path = &summary.populations.local_path;
    let populations_file_name = populations_path.file_name().unwrap().to_str().unwrap();
    let alignment_reader = fasta::Reader::from_file(populations_path)
        .map_err(|e| eyre!(e))
        .wrap_err("Failed to read file: {populations_path:?}")?;

    // keep track of population names in alignment, cross-reference against
    // lineage notes + alias_key later
    let mut alignment_populations = Vec::new();
    for result in alignment_reader.records() {
        let record =
            result.wrap_err(eyre!("Failed to parse file: {populations_path:?}"))?;
        let lineage = record.id().to_string();
        alignment_populations.push(lineage);
    }

    // ------------------------------------------------------------------------
    // Parent Child Relationships
    // ------------------------------------------------------------------------

    let mut graph_order = Vec::new();
    let mut graph_data = BTreeMap::new();

    for row in lineage_table.rows {
        let lineage = row[lineage_col_i].to_string();

        // Lineages that start with '*' have been withdrawn
        if lineage.starts_with('*') || lineage == String::new() {
            continue;
        }

        let parents = get_lineage_parents(&lineage, &alias_key)?;
        graph_order.push(lineage.clone());
        graph_data.insert(lineage, parents);
    }

    // ------------------------------------------------------------------------
    // Graph (Phylogeny)
    // ------------------------------------------------------------------------

    let mut phylogeny = Phylogeny::new();

    // Add root node
    let name = "root".to_string();
    phylogeny.graph.add_node(name);

    // todo!() Do this twice? in case lineages are accidentally out of order?

    // Add descendants
    for name in graph_order {
        let id = phylogeny.graph.add_node(name.clone());
        if !graph_data.contains_key(&name) {
            return Err(
                eyre!("Parents of {name} are unknown in the phylogeny graph.")
                    .suggestion(
                        "Could the lineage_notes be out of order or misformatted?",
                    )
                    .suggestion("Parents are required to appear before children."),
            );
        }
        let parents = graph_data.get(&name).unwrap();

        debug!("Population: {name}; Parents: {parents:?}");

        // If multiple parents add this to recombinants list
        if parents.len() > 1 {
            phylogeny.recombinants.push(name.clone());
        }
        for parent in parents {
            if phylogeny.get_node(parent).is_err() {
                return Err(eyre!("Parental lineage {parent} is not in the graph.")
                    .suggestion(
                        "Are the alias_key.json and lineage_notes.txt out of sync?",
                    )
                    .suggestion("Please check if {parent} is in the alias key."));
            }
            let parent_id = phylogeny.get_node(parent)?;
            phylogeny.graph.add_edge(parent_id, id, 1);
        }
    }

    phylogeny.recombinants_all = phylogeny.get_recombinants_all()?;

    // ------------------------------------------------------------------------
    // Consistency Check
    // ------------------------------------------------------------------------

    // All population names in all files
    let mut all_populations = alignment_populations.clone();
    all_populations.extend(notes_lineages.clone());
    all_populations.extend(alias_key.keys().cloned().collect_vec());
    all_populations.sort();
    all_populations.dedup();

    // create table to store consistency info
    let inconsistency_table_path = output_dir.join("inconsistency.tsv");
    let mut inconsistency_table = Table::new();

    info!("Writing dataset inconsistency table: {inconsistency_table_path:?}");

    inconsistency_table.headers = vec!["population", "present", "absent"]
        .into_iter()
        .map(String::from)
        .collect_vec();

    let population_col_i = inconsistency_table.header_position("population")?;
    let present_col_i = inconsistency_table.header_position("present")?;
    let absent_col_i = inconsistency_table.header_position("absent")?;

    // missing population sequences
    let mut missing_pop_seq = Vec::new();

    // check consistency between populations fasta, lineage_notes, alias_key, phylogeny
    all_populations.iter().for_each(|pop| {
        let mut present_file_names = Vec::new();
        let mut absent_file_names = Vec::new();

        // Alias Key
        if alias_key.contains_key(pop) {
            present_file_names.push(&alias_key_file_name);
        } else if !pop.contains('.') {
            absent_file_names.push(&alias_key_file_name);
        }

        // Lineage Notes
        if notes_lineages.contains(pop) {
            present_file_names.push(&lineage_notes_file_name);
        } else {
            absent_file_names.push(&lineage_notes_file_name);
        }

        // Populations Fasta
        if alignment_populations.contains(pop) {
            present_file_names.push(&populations_file_name);
        } else {
            absent_file_names.push(&populations_file_name);
            missing_pop_seq.push(pop.clone());
        }

        // Phylogeny
        // self.graph.node_references()

        if !absent_file_names.is_empty() {
            let mut row = vec![String::new(); inconsistency_table.headers.len()];
            row[population_col_i] = pop.to_string();
            row[present_col_i] = present_file_names.iter().join(",");
            row[absent_col_i] = absent_file_names.iter().join(",");
            inconsistency_table.rows.push(row);
        }
    });

    inconsistency_table.write(&inconsistency_table_path)?;

    Ok(phylogeny)
}

/// Get all immediate parents of a SARS-CoV-2 lineage.
pub fn get_lineage_parents(
    lineage: &str,
    alias_key: &BTreeMap<String, Vec<String>>,
) -> Result<Vec<String>, Report> {
    let mut parents: Vec<String> = Vec::new();

    // If Recombinant with multiple parents, if so, it will be in the alias
    // key with the parents listed.
    if alias_key.contains_key(lineage) {
        let lineage_paths = &alias_key[lineage];
        if lineage_paths.len() > 1 {
            // Dedup in case multiple breakpoints/parents
            parents = lineage_paths.clone().into_iter().unique().collect();
            return Ok(parents);
        }
    }

    // Otherwise, single parent
    let decompress = decompress_lineage(lineage, alias_key).unwrap();

    // Ex. BA.5.2 -> ["BA", "5", "2"]
    let decompress_parts =
        decompress.split('.').map(|p| p.to_string()).collect::<Vec<_>>();

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

/// Compress a SARS-CoV-2 lineage name into it's aliased form.
pub fn compress_lineage(
    lineage: &String,
    alias_key: &BTreeMap<String, Vec<String>>,
) -> Result<String, Report> {
    // By default, set compression level to self
    let mut compress = lineage.to_string();

    // Reverse the alias-> lineage path lookup
    let mut alias_key_rev: BTreeMap<String, String> = BTreeMap::new();

    for (alias, lineage_paths) in alias_key {
        // Skip over recombinants with multiple parents, don't need their lookup
        if lineage_paths.len() > 1 {
            continue;
        }
        let lineage_path = lineage_paths[0].clone();
        alias_key_rev.insert(lineage_path, alias.clone());
    }

    // Ex. BA.5.2 -> ["BA", "5", "2"]
    let compress_parts = compress.split('.').map(|p| p.to_string()).collect::<Vec<_>>();

    if compress_parts.len() > 1 {
        for i in (0..compress_parts.len()).rev() {
            let compress_subset = compress_parts[0..i].join(".");

            if alias_key_rev.contains_key(&compress_subset) {
                compress = alias_key_rev.get(&compress_subset).unwrap().clone();
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

/// Decompress a SARS-CoV-2 lineage name into it's full unaliased form.
///
/// # Arguments
///
///  * `lineage` | A string slice that contains a SARS-CoV-2 lineage name.
///  * `alias_key` | A mapping of SARS-CoV-2 aliases to decompressed lineage paths.
/// ```
pub fn decompress_lineage(
    lineage: &str,
    alias_key: &BTreeMap<String, Vec<String>>,
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
    let lineage_paths = alias_key.get(&lineage_prefix).cloned().unwrap_or_default();
    // Not multiple recombinant parents
    if lineage_paths.len() == 1 {
        decompress = lineage_paths[0].clone();
        // Add back our suffix numbers
        if lineage_level > 1 {
            decompress = format!("{decompress}.{lineage_suffix}");
        }
    }

    Ok(decompress)
}

/// Read SARS-CoV-2 alias key from file.
///
/// Reads the JSON file into a Map where keys are lineage names, and
/// values are a vector of parents.
pub fn read_alias_key(path: &Path) -> Result<BTreeMap<String, Vec<String>>, Report> {
    // read to json
    let alias_key_str =
        std::fs::read_to_string(path).expect("Couldn't read alias_key file.");
    let alias_key_val: serde_json::Value =
        serde_json::from_str(&alias_key_str).expect("Couldn't convert alias_key to json");
    // deserialize to object (raw, mixed types)
    let alias_key_raw: serde_json::Map<String, serde_json::Value> = alias_key_val
        .as_object()
        .expect("Couldn't convert alias_key json to json Map")
        .clone();

    // This should probably be a custom serializer fn for brevity,
    //   but I don't know how to do that yet :)
    let mut alias_key: BTreeMap<String, Vec<String>> = BTreeMap::new();

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
