pub mod remote_file;
pub mod table;

use crate::dataset::attributes::Tag;
use crate::utils::remote_file::RemoteFile;
use chrono::prelude::*;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use color_eyre::Help;
use itertools::Itertools;
use log::{debug, warn};
use reqwest::header::{ACCESS_CONTROL_EXPOSE_HEADERS, USER_AGENT};
use std::collections::BTreeMap;
use std::fs::{remove_file, write, File};
use std::io::Read;
use std::path::{Path, PathBuf};
use tempfile::TempDir;
use zstd::stream::read::Decoder;

/// Download file from url to path, with optional decompression.
pub async fn download_file(
    url: &str,
    output_path: &Path,
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

pub fn check_github_response(response: &reqwest::Response) -> Result<(), Report> {
    let url = response.url().to_string();

    if !response.status().is_success() {
        // --------------------------------------------------------------------
        // Check if the API rate limit was exceeded

        // todo!() this is some pretty risky unwrapping here
        let rate_limit_remaining: u32 = response
            .headers()
            .get("x-ratelimit-remaining")
            .unwrap()
            .to_str()?
            .parse()
            .unwrap();
        if rate_limit_remaining == 0 {
            let rate_limit_reset: i64 = response
                .headers()
                .get("x-ratelimit-reset")
                .unwrap()
                .to_str()?
                .parse()
                .unwrap();
            let rate_limit_reset: DateTime<Local> =
                DateTime::<Utc>::from_timestamp(rate_limit_reset, 0)
                    .expect("invalid timestamp")
                    .into();

            return Err(
                eyre!("GitHub API rate limit has been exceeded.")
                .suggestion(format!("Please wait for the rate limit to reset at: {rate_limit_reset:?}"))
                .suggestion("Alternatively, set the environment variables GITHUB_USERNAME and GITHUB_TOKEN.")
                .suggestion("https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens")
            );
        }
        // --------------------------------------------------------------------
        // Other unhandled errors
        else {
            return Err(eyre!(
                "GitHub query had status code {status}: {url}",
                status = response.status()
            ));
        }
    }

    Ok(())
}

/// Query and download files using the GitHub API
pub async fn download_github(
    repo: &str,
    tag: &Tag,
    remote_path: &str,
    output_path: &Path,
    sha: &Option<String>,
) -> Result<RemoteFile, Report> {
    // GitHub API Setup
    let github_token: Option<String> = std::env::var("GITHUB_TOKEN").ok();
    let github_username = std::env::var("GITHUB_USERNAME").unwrap_or("".to_string());
    let github_api_version = "2022-11-28";
    let user_agent = format!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
    let client = reqwest::Client::new();

    // GitHub API Query
    let mut query = vec![("path", remote_path), ("per_page", "1"), ("page", "1")];
    // if manual SHA was specified
    if let Some(sha) = sha {
        query.push(("sha", sha));
    }

    // convert to string
    let mut query =
        query.into_iter().map(|(k, v)| (String::from(k), String::from(v))).collect_vec();

    // Note: the reqwest::RequestBuilder doesn't implement Clone for the
    // non-blocking (asynchronous) version :( we're going to have to full
    // define the request over and over again.

    // --------------------------------------------------------------------------
    // STEP 1: Archive Pagination

    // Use the Archive Tag as a maximum date filter (&until=...)
    if matches!(tag, Tag::Archive(_)) {
        query.push(("until".to_string(), tag.to_string()));
    }

    let request = client
        .get(format!("https://api.github.com/repos/{repo}/commits"))
        .query(&query)
        .header(USER_AGENT, &user_agent)
        .header(ACCESS_CONTROL_EXPOSE_HEADERS, "Link")
        .header("X-GitHub-Api-Version", github_api_version)
        .basic_auth(&github_username, github_token.clone());
    let response = request.send().await?;
    check_github_response(&response)?;

    let url = response.url().to_string();

    // extract the "sha" and "date" key from the json body
    let body: Vec<BTreeMap<String, serde_json::Value>> = response.json().await?;
    if body.is_empty() {
        return Err(eyre!("No GitHub commits were found for: {}", url)
            .suggestion(format!(
                "Does your dataset tag ({tag}) predate the creation of this file?"
            ))
            .suggestion(format!(
                "Repository: https://github.com/{repo}, File: {remote_path:?}"
            )));
    }

    let sha = body[0]["sha"].to_string().replace('"', "");
    let commit_date = body[0]["commit"]["author"]["date"].to_string().replace('"', "");
    let date_created: DateTime<Utc> = DateTime::parse_from_rfc3339(&commit_date)?.into();

    // --------------------------------------------------------------------------
    // STEP 2: DOWNLOAD

    let download_url =
        format!("https://raw.githubusercontent.com/{repo}/{sha}/{remote_path}");

    // Identify decompression mode
    // TBD! todo!() make this an enum of implemented decompression types
    let ext = path_to_ext(Path::new(&download_url))?;
    let decompress = ext == "zst";

    // Download the file
    debug!("Downloading file: {download_url} to {output_path:?}");
    download_file(&download_url, output_path, decompress).await?;

    // Store all the information about the remote file for the dataset summary
    let remote_file = RemoteFile {
        url: download_url,
        sha,
        local_path: output_path.to_path_buf(),
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
    let delim = match ext.as_str() {
        "tsv" => '\t',
        "csv" => ',',
        "txt" => {
            warn!("File extension .txt is assumed to be tab-delimited: {path:?}");
            '\t'
        }
        _ => {
            return Err(eyre!("Unknown file extension: {ext:?}")
                .suggestion("Options are tsv or csv."))
        }
    };

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
