use crate::dataset::attributes::Annotations;
use crate::dataset::edge_cases::EdgeCase;
use crate::utils;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use itertools::Itertools;
use log::{debug, info};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

// ----------------------------------------------------------------------------
// URLs
// ----------------------------------------------------------------------------

pub const POPULATIONS_URL: &str = "https://raw.githubusercontent.com/corneliusroemer/pango-sequences/main/data/pango-consensus-sequences_genome-nuc.fasta.zst";
pub const REFERENCE_URL: &str = "https://raw.githubusercontent.com/nextstrain/ncov/master/data/references_sequences.fasta";
pub const ALIAS_KEY_URL: &str = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json";
pub const LINEAGE_NOTES_URL: &str = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt";

// ----------------------------------------------------------------------------
// Downloaders
// ----------------------------------------------------------------------------

/// Download the SARS-CoV-2 lineage notes.
///
/// The lineage notes has two columns: 'Lineage', 'Description'.
/// We only need the 'Lineage' column, to get the full list of all lineages.
pub async fn download_lineage_notes(dataset_dir: &Path) -> Result<PathBuf, Report> {
    let decompress = false;
    let output_path = dataset_dir.join("lineage_notes.txt");
    utils::download_file(LINEAGE_NOTES_URL, &output_path, decompress)
        .await
        .wrap_err_with(|| eyre!("Unable to download alias lineage notes."))?;
    Ok(output_path)
}

/// Download the SARS-CoV-2 alias key.
///
/// The alias key is a JSON mapping lineage names to their parents.
/// Needed to construct the phylogeny and identify known recombinants.
pub async fn download_alias_key(dataset_dir: &Path) -> Result<PathBuf, Report> {
    let decompress = false;
    let output_path = dataset_dir.join("alias_key.json");
    utils::download_file(ALIAS_KEY_URL, &output_path, decompress)
        .await
        .wrap_err_with(|| eyre!("Unable to download alias from key."))?;

    Ok(output_path)
}

// ----------------------------------------------------------------------------
// Loaders
// ----------------------------------------------------------------------------

/// Load the SARS-CoV-2 alias key.
///
/// Reads the JSON file into a Map where keys are lineage names, and
/// values are a vector of parents.
pub fn load_alias_key(path: &Path) -> Result<BTreeMap<String, Vec<String>>, Report> {
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

// ----------------------------------------------------------------------------
// Parsers
// ----------------------------------------------------------------------------

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

/// Decompress a SARS-CoV-2 lineage name into it's full unaliased form.
///
/// # Arguments
///
///  * `lineage` | `String` | A String of a SARS-CoV-2 lineage name.
///  * `alias_key` | `&BTreeMap<String, Vec<String>>` | A mapping of SARS-CoV-2 aliases to decompressed lineage paths.
///
/// # Basic usage:
///
/// ```
/// decompress_lineage("BA.5.2", alias_key)?;
/// ```
///
/// # Example
///
/// ```no run
/// let alias_key = import_alias_key(".")?;
/// decompress_lineage("BA.5.2", alias_key)?;
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

// ----------------------------------------------------------------------------
// Creators
// ----------------------------------------------------------------------------

/// Create SARS-CoV-2 genome annotations.
pub fn create_annotations() -> Result<Annotations, Report> {
    let annotations = Annotations {
        gene: vec![
            "ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8",
            "ORF9b",
        ]
        .into_iter()
        .map(String::from)
        .collect_vec(),
        abbreviation: vec!["1a", "1b", "S", "3a", "E", "M", "6", "7a", "7b", "8", "9b"]
            .into_iter()
            .map(String::from)
            .collect_vec(),
        start: vec![
            266, 13468, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28284,
        ],
        end: vec![
            13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 28577,
        ],
    };

    Ok(annotations)
}

/// Create SARS-CoV-2 recombinant edge cases.
pub fn create_edge_cases() -> Result<Vec<EdgeCase>, Report> {
    info!("Creating SARS-CoV-2 edge cases.");
    let mut edge_cases: Vec<EdgeCase> = Vec::new();

    // --------------------------------------------------------------------
    // XCF
    // XCF is XBB and FE.1 (XBB.1.18.1) with no unique subs from XBB

    debug!("Creating edge case: XCF");
    let xcf = EdgeCase {
        population: "XCF".to_string(),
        max_parents: 2,
        max_iter: 1,
        min_consecutive: 1,
        min_length: 1,
        min_subs: 0,
        unbiased: false,
    };
    edge_cases.push(xcf);

    // --------------------------------------------------------------------
    // Finish

    Ok(edge_cases)
}

/// Create SARS-CoV-2 phylogeny graph data.
pub async fn create_graph_data(
    dataset_dir: &Path,
) -> Result<(BTreeMap<String, Vec<String>>, Vec<String>), Report> {
    // ------------------------------------------------------------------------
    // Download Data

    // Download the alias key
    let alias_key_path = download_alias_key(dataset_dir).await?;
    // Download the lineage notes
    let lineage_notes_path = download_lineage_notes(dataset_dir).await?;

    // ------------------------------------------------------------------------
    // Load Data

    // read alias key into Map
    let alias_key = load_alias_key(&alias_key_path)?;
    // read lineage notes into Table
    let lineage_table = utils::read_table(&lineage_notes_path)?;
    // identify which column is 'Lineage'
    let lineage_col_i = lineage_table.header_position("Lineage")?;

    // ------------------------------------------------------------------------
    // Map parent-child relationships

    let mut graph_data: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut graph_order: Vec<String> = Vec::new();

    for row in lineage_table.rows {
        let lineage = row[lineage_col_i].to_string();

        // Lineages that start with '*' have been withdrawn
        if lineage.starts_with('*') {
            continue;
        }

        let parents = get_lineage_parents(&lineage, &alias_key)?;

        graph_order.push(lineage.clone());
        graph_data.insert(lineage, parents);
    }

    Ok((graph_data, graph_order))
}