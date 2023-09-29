use crate::cli;
use crate::dataset::edge_cases::EdgeCase;
use crate::utils;
use crate::utils::{remote_file::RemoteFile, table::Table};
use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use itertools::Itertools;
use log::{debug, warn};
use std::collections::BTreeMap;
use std::path::Path;

// ----------------------------------------------------------------------------
// Downloaders
// ----------------------------------------------------------------------------

pub async fn download_reference(
    args: &cli::DatasetDownloadArgs,
) -> Result<RemoteFile, Report> {
    let repo = "nextstrain/ncov";
    let remote_path = "data/references_sequences.fasta";
    let output_path = args.output_dir.join("reference.fasta");
    let sha: Option<String> = None;
    let remote_file =
        utils::download_github(repo, &args.tag, remote_path, &output_path, &sha)
            .await
            .wrap_err_with(|| eyre!("Unable to download sars-cov-2 reference fasta."))?;
    Ok(remote_file)
}
pub async fn download_populations(
    args: &cli::DatasetDownloadArgs,
) -> Result<RemoteFile, Report> {
    let repo = "corneliusroemer/pango-sequences";
    let remote_path = "data/pango-consensus-sequences_genome-nuc.fasta.zst";
    let output_path = args.output_dir.join("populations.fasta");
    let sha: Option<String> = None;
    let remote_file =
        utils::download_github(repo, &args.tag, remote_path, &output_path, &sha)
            .await
            .wrap_err_with(|| {
                eyre!("Unable to download sars-cov-2 populations fasta.")
            })?;
    Ok(remote_file)
}

/// Download the SARS-CoV-2 alias key.
///
/// The alias key is a JSON mapping lineage names to their parents.
/// Needed to construct the phylogeny and identify known recombinants.
pub async fn download_alias_key(
    args: &cli::DatasetDownloadArgs,
) -> Result<RemoteFile, Report> {
    let repo = "cov-lineages/pango-designation";
    let remote_path = "pango_designation/alias_key.json";
    let output_path = args.output_dir.join("alias_key.json");
    let sha: Option<String> = None;
    let remote_file =
        utils::download_github(repo, &args.tag, remote_path, &output_path, &sha)
            .await
            .wrap_err_with(|| eyre!("Unable to download sars-cov-2 alias key."))?;
    Ok(remote_file)
}

/// Download the SARS-CoV-2 lineage notes.
///
/// The lineage notes has two columns: 'Lineage', 'Description'.
/// We only need the 'Lineage' column, to get the full list of all lineages.
pub async fn download_lineage_notes(
    args: &cli::DatasetDownloadArgs,
) -> Result<RemoteFile, Report> {
    let repo = "cov-lineages/pango-designation";
    let remote_path = "lineage_notes.txt";
    let output_path = args.output_dir.join("lineage_notes.txt");
    let sha: Option<String> = None;
    let remote_file =
        utils::download_github(repo, &args.tag, remote_path, &output_path, &sha)
            .await
            .wrap_err_with(|| eyre!("Unable to download sars-cov-2 lineage notes."))?;
    Ok(remote_file)
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
pub fn create_annotations() -> Result<Table, Report> {
    let mut table = Table::new();

    let headers = vec!["gene", "abbreviation", "start", "end"];
    let rows = vec![
        vec!["ORF1a", "1a", "266", "13468"],
        vec!["ORF1b", "1b", "13468", "21555"],
        vec!["S", "S", "21563", "25384"],
        vec!["ORF3a", "3a", "25393", "26220"],
        vec!["E", "E", "26245", "26472"],
        vec!["M", "M", "26523", "27191"],
        vec!["ORF6", "6", "27202", "27387"],
        vec!["ORF7a", "7a", "27394", "27759"],
        vec!["ORF7b", "7b", "27756", "27887"],
        vec!["ORF8", "8", "27894", "28259"],
        vec!["ORF9b", "9b", "28284", "28577"],
    ];

    // Convert values to String
    table.headers = headers.into_iter().map(String::from).collect_vec();
    table.rows = rows
        .into_iter()
        .map(|row| row.into_iter().map(String::from).collect_vec())
        .collect_vec();

    Ok(table)
}

/// Create SARS-CoV-2 recombinant edge cases.
pub fn create_edge_cases() -> Result<Vec<EdgeCase>, Report> {
    let mut edge_cases: Vec<EdgeCase> = Vec::new();

    // --------------------------------------------------------------------
    // B
    // B is potentially the reference itself, with no mutations?

    debug!("Creating edge case: B");
    let mut b = EdgeCase::new();
    b.population = "B".to_string();

    // --------------------------------------------------------------------
    // XCF
    // XCF is XBB and FE.1 (XBB.1.18.1) with no unique subs from XBB

    debug!("Creating edge case: XCF");
    let xcf = EdgeCase {
        population: "XCF".to_string(),
        max_parents: 2,
        max_iter: 1,
        min_consecutive: 2,
        min_length: 500,
        min_subs: 0,
        unbiased: false,
    };
    edge_cases.push(xcf);

    // --------------------------------------------------------------------
    // XCG
    // XCG is BA.5.2 and XBB.1 with only 2 consecutive bases

    debug!("Creating edge case: XCG");
    let xcg = EdgeCase {
        population: "XCG".to_string(),
        max_parents: 2,
        max_iter: 1,
        min_consecutive: 2,
        min_length: 500,
        min_subs: 1,
        unbiased: false,
    };
    edge_cases.push(xcg);

    // --------------------------------------------------------------------
    // Finish

    Ok(edge_cases)
}

/// Create SARS-CoV-2 phylogeny graph data.
pub async fn create_graph_data(
    args: &cli::DatasetDownloadArgs,
) -> Result<(BTreeMap<String, Vec<String>>, Vec<String>), Report> {
    // ------------------------------------------------------------------------
    // Load Data

    // read alias key into Map
    let alias_key_path = args.output_dir.join("alias_key.json");
    let alias_key = load_alias_key(&alias_key_path)?;

    // read lineage notes into Table
    let lineage_notes_path = args.output_dir.join("lineage_notes.txt");
    let lineage_table = utils::read_table(&lineage_notes_path)?;
    // identify which column is 'Lineage'
    let lineage_col_i = lineage_table.header_position("Lineage")?;
    // get a list of lineages in the notes
    let notes_lineages = lineage_table
        .rows
        .iter()
        .filter(|row| !row[lineage_col_i].starts_with('*'))
        .map(|row| row[lineage_col_i].to_string())
        .collect_vec();

    // read populations fasta, to check if any lineages are missing in notes
    let populations_path = args.output_dir.join("populations.fasta");
    let alignment_reader = fasta::Reader::from_file(&populations_path)
        .map_err(|e| eyre!(e))
        .wrap_err("Unable to read populations: {populations_path:?}")?;

    for result in alignment_reader.records() {
        let record =
            result.wrap_err(eyre!("Unable to parse alignment: {populations_path:?}"))?;
        let lineage = record.id().to_string();
        if !notes_lineages.contains(&lineage) {
            warn!("Lineage {lineage} is in {populations_path:?} but not in {lineage_notes_path:?}.");
        }
    }

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
