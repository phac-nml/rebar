pub mod table;

use crate::utils::table::Table;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use color_eyre::Help;
use itertools::Itertools;
use log::warn;
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
