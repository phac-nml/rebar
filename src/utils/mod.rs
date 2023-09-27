pub mod remote_file;
pub mod table;

use crate::dataset::attributes::Tag;
use crate::utils::remote_file::RemoteFile;
use crate::utils::table::Table;
use chrono::prelude::*;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use color_eyre::Help;
use itertools::Itertools;
use log::{debug, info, warn};
use reqwest::header::{ACCESS_CONTROL_EXPOSE_HEADERS, USER_AGENT};
use std::collections::BTreeMap;
use std::fs::{remove_file, write, File};
use std::io::{self, BufRead, Read};
use std::path::{Path, PathBuf};
use tempfile::TempDir;
use zstd::stream::read::Decoder;

/// Download file from url to path, with optional decompression.
pub async fn download_file(
    url: &str,
    output_path: &PathBuf,
    decompress: bool,
) -> Result<(), Report> {
    let ext = Path::new(&url).extension().unwrap().to_str().unwrap();

    let response = reqwest::get(url).await?;
    if response.status() != 200 {
        return Err(eyre!(
            "Unable to download file: {url}\nStatus code {}.",
            response.status()
        ));
    }

    if decompress {
        // Write bytes to a tmp file
        let tmp_dir = TempDir::new()?;
        let tmp_path = PathBuf::from(tmp_dir.path()).join(format!("tmpfile.{ext}"));
        let content = response.bytes().await?;
        write(&tmp_path, content)
            .wrap_err_with(|| eyre!("Unable to write file: {tmp_path:?}"))?;
        decompress_file(&tmp_path, output_path, true)?;
    } else {
        let content = response.text().await?;
        write(output_path, content)
            .wrap_err_with(|| eyre!("Unable to write file: {output_path:?}"))?;
    }

    Ok(())
}

/// Query and download files using the GitHub API
pub async fn download_github(
    repo: &str,
    tag: &Tag,
    remote_path: &str,
    output_path: &PathBuf,
) -> Result<RemoteFile, Report> {
    let query_url = format!(
        "https://api.github.com/repos/{repo}/commits?per_page=1&path={remote_path}"
    );
    let download_url = format!(
        "https://raw.githubusercontent.com/{repo}/{sha}/{remote_path}",
        sha = "{sha}"
    );

    let client = reqwest::Client::new();
    let user_agent = format!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));

    // --------------------------------------------------------------------------
    // STEP 1: Commit SHA

    let query_url = match tag {
        // For Unknown tag, error out
        Tag::Unknown => {
            return Err(eyre!("GitHub download for tag {tag:?} is not implemented."))
        }
        // For the latest tag, just use query_url plainly
        Tag::Latest => query_url,
        // For archive, we need to deal with query pagination of historical commits
        Tag::Archive(_) => {
            let tag_date: DateTime<Utc> =
                DateTime::parse_from_rfc3339(&tag.to_string())?.into();
            let mut min_date_query_url =
                format!("{query_url}&since={}", tag_date.format("%Y-%m-%d"));
            let response = client
                .get(&min_date_query_url)
                .header(USER_AGENT, &user_agent)
                .header(ACCESS_CONTROL_EXPOSE_HEADERS, "Link")
                .send()
                .await?;
            let headers = response.headers();

            // extract the "last" page url, to get the oldest commit on or after tag date
            // link might not be present if no pagination was needed (only 1 result)
            if headers.contains_key("link") {
                let link = headers["link"].to_str()?;
                // link format:
                //  [0] "<https://api.github.com/repositories/538600532/commits?since=2023-08-17&per_page=1&page=2>;"
                //  [1]: rel=\"next\",
                //  [2] <https://api.github.com/repositories/538600532/commits?since=2023-08-17&per_page=1&page=99>;
                //  [3] <rel=\"last\""
                let link_parts = link.split(' ').map(|s| s.to_string()).collect_vec();
                min_date_query_url = link_parts[2].replace(['<', '>', ';', ' '], "");
                min_date_query_url
            } else {
                query_url
            }
        }
    };

    let response = client
        .get(&query_url)
        .header(USER_AGENT, &user_agent)
        .send()
        .await?;

    // extract the "sha" and "date" key from the json body
    let body: Vec<BTreeMap<String, serde_json::Value>> = response.json().await?;
    if body.is_empty() {
        return Err(
            eyre!("No GitHub commits were found for the following query.")
                .suggestion(format!("GitHub API query: {query_url}")),
        );
    }

    let sha = body[0]["sha"].to_string().replace('"', "");
    let commit_date = body[0]["commit"]["author"]["date"]
        .to_string()
        .replace('"', "");
    let date_created: DateTime<Utc> = DateTime::parse_from_rfc3339(&commit_date)?.into();

    // --------------------------------------------------------------------------
    // STEP 2: DOWNLOAD

    let download_url = download_url.replace("{sha}", &sha);

    // Identify decompression mode
    let ext = path_to_ext(Path::new(&download_url))?;
    let decompress = ext != "fasta" && ext != "fa";

    // Download the file
    info!("Downloading file: {download_url} to {output_path:?}");
    download_file(&download_url, output_path, decompress).await?;

    // Store all the information about the remote file for the dataset summary
    let remote_file = RemoteFile {
        url: download_url,
        sha,
        local_path: output_path.clone(),
        date_created,
        date_downloaded: Utc::now(),
    };
    debug!("Downloaded file: {remote_file:?}");

    Ok(remote_file)
}

/// Decompress file, optionally inplace
pub fn decompress_file(input: &Path, output: &Path, inplace: bool) -> Result<(), Report> {
    let ext = input
        .extension()
        .ok_or_else(|| eyre!("Unable to parse extension from file: {input:?}"))?;

    match ext.to_str().unwrap() {
        "zst" => {
            let reader = File::open(input)?;
            let mut decoder = Decoder::new(reader)?;
            let mut buffer = String::new();
            decoder.read_to_string(&mut buffer)?;
            write(output, buffer)
                .wrap_err(format!("Unable to write file: {:?}", output))?;

            if inplace {
                remove_file(input)?;
            }
        }
        _ => return Err(eyre!("Decompression for .{ext:?} is not implemented yet.")),
    };

    Ok(())
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
/// Source: https://doc.rust-lang.org/rust-by-example/std_misc/file/read_lines.html
pub fn read_lines(path: &Path) -> Result<io::Lines<io::BufReader<File>>, Report> {
    // attempt to open the file path
    let file = File::open(path)
        .wrap_err_with(|| format!("Failed to open table at: {path:?}"))?;
    // read in the lines
    let lines = io::BufReader::new(file).lines();
    Ok(lines)
}

pub fn read_table(path: &Path) -> Result<Table, Report> {
    let mut table = Table::new();

    // lookup delimiter from file extension
    let delim = path_to_delim(path)?;

    for line in (read_lines(path)?).flatten() {
        let row = line
            .split(delim)
            .collect_vec()
            .into_iter()
            .map(String::from)
            .collect_vec();
        // if headers are empty, this is the first line, write headers
        if table.headers.is_empty() {
            table.headers = row;
        }
        // otherwise regular row
        else {
            table.rows.push(row);
        }
    }

    table.path = path.to_path_buf();

    Ok(table)
}

pub fn ext_to_delim(ext: &str) -> Result<char, Report> {
    let delim = match ext {
        "tsv" => '\t',
        "csv" => ',',
        "txt" => {
            warn!("File extension .txt is assumed to be tab-delimited.");
            '\t'
        }
        _ => {
            return Err(eyre!("Unknown file extension: {ext:?}")
                .suggestion("Options are tsv or csv."))
        }
    };

    Ok(delim)
}

pub fn path_to_delim(path: &Path) -> Result<char, Report> {
    // get the path extension
    let ext = path_to_ext(path)?;

    // convert extension to the expected delimiter
    let delim = ext_to_delim(&ext)?;

    Ok(delim)
}

pub fn path_to_ext(path: &Path) -> Result<String, Report> {
    let result = path.extension();
    let ext = match result {
        Some(ext) => ext.to_os_string().into_string().unwrap(),
        None => return Err(eyre!("Unable to parse extension from file: {path:?}")),
    };

    Ok(ext)
}
