use crate::dataset::attributes::Tag;
use crate::utils::{download_github, remote_file::RemoteFile};
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use std::path::Path;

pub async fn reference(tag: &Tag, output_path: &Path) -> Result<RemoteFile, Report> {
    let repo = "nextstrain/ncov";
    let remote_path = "data/references_sequences.fasta";
    let sha: Option<String> = None;
    let remote_file = download_github(repo, tag, remote_path, output_path, &sha)
        .await
        .wrap_err_with(|| eyre!("Failed downloading sars-cov-2 reference fasta."))?;
    Ok(remote_file)
}
pub async fn populations(tag: &Tag, output_path: &Path) -> Result<RemoteFile, Report> {
    let repo = "corneliusroemer/pango-sequences";
    let remote_path = "data/pango-consensus-sequences_genome-nuc.fasta.zst";
    let sha: Option<String> = None;
    let remote_file = download_github(repo, tag, remote_path, output_path, &sha)
        .await
        .wrap_err_with(|| eyre!("Failed downloading sars-cov-2 populations fasta."))?;
    Ok(remote_file)
}

/// Download the SARS-CoV-2 alias key.
///
/// The alias key is a JSON mapping lineage names to their parents.
/// Needed to construct the phylogeny and identify known recombinants.
pub async fn alias_key(tag: &Tag, output_path: &Path) -> Result<RemoteFile, Report> {
    let repo = "cov-lineages/pango-designation";
    let remote_path = "pango_designation/alias_key.json";
    let sha: Option<String> = None;
    let remote_file = download_github(repo, tag, remote_path, output_path, &sha)
        .await
        .wrap_err_with(|| eyre!("Failed downloading sars-cov-2 alias key."))?;
    Ok(remote_file)
}

/// Download the SARS-CoV-2 lineage notes.
///
/// The lineage notes has two columns: 'Lineage', 'Description'.
/// We only need the 'Lineage' column, to get the full list of all lineages.
pub async fn lineage_notes(tag: &Tag, output_path: &Path) -> Result<RemoteFile, Report> {
    let repo = "cov-lineages/pango-designation";
    let remote_path = "lineage_notes.txt";
    let sha: Option<String> = None;
    let remote_file = download_github(repo, tag, remote_path, output_path, &sha)
        .await
        .wrap_err_with(|| eyre!("Failed downloading sars-cov-2 lineage notes."))?;
    Ok(remote_file)
}

/// Download the SARS-CoV-2 nameTable mapping clades to lineage names.
pub async fn clade_to_lineage(
    tag: &Tag,
    output_path: &Path,
) -> Result<RemoteFile, Report> {
    let repo = "hodcroftlab/covariants";
    let remote_path = "web/data/nameTable.json";
    let sha: Option<String> = None;
    let remote_file =
        download_github(repo, tag, remote_path, output_path, &sha).await.wrap_err_with(
            || eyre!("Failed downloading sars-cov-2 clade_to_lineage nameTable."),
        )?;
    Ok(remote_file)
}
