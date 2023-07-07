use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use csv;
use itertools::Itertools;
use std::fs::{remove_file, write, File};
use std::io::Read;
use std::path::{Path, PathBuf};
use tempfile::TempDir;
use zstd::stream::read::Decoder;

#[derive(Debug, Clone)]
pub struct Table {
    pub headers: Vec<String>,
    pub rows: Vec<Vec<String>>,
    pub path: PathBuf,
}

impl Default for Table {
    fn default() -> Self {
        Self::new()
    }
}

impl Table {
    pub fn new() -> Self {
        Table {
            path: PathBuf::new(),
            headers: Vec::new(),
            rows: Vec::new(),
        }
    }

    pub fn from_tsv(path: &Path) -> Result<Self, Report> {
        let mut table = Table::new();

        // init reader from file path
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .wrap_err_with(|| format!("Failed to open table at: {}", path.display()))?;

        // read in headers and convert to String
        let headers = reader.headers()?.iter().map(String::from).collect_vec();

        // read in records, then convert to rows of String
        let records: Vec<csv::StringRecord> = reader
            .records()
            .collect::<Result<Vec<csv::StringRecord>, csv::Error>>()?;
        let rows: Vec<Vec<String>> = records
            .iter()
            .map(|row| row.iter().map(String::from).collect_vec())
            .collect_vec();

        table.headers = headers;
        table.rows = rows;
        table.path = path.to_path_buf();

        Ok(table)
    }

    pub fn header_position(&self, header: &str) -> Result<usize, Report> {
        let pos = self
            .headers
            .iter()
            .position(|h| h == header)
            .ok_or_else(|| {
                eyre!("Column '{header}' was not found in table: {:?}.", self.path)
            })?;

        Ok(pos)
    }

    /// write to file
    pub fn write(&self, output_path: &Path, delim: Option<char>) -> Result<(), Report> {
        let delim = match delim {
            Some(c) => c as u8,
            None => b'\t',
        };
        let mut writer = csv::WriterBuilder::new()
            .delimiter(delim)
            .from_path(output_path)
            .wrap_err_with(|| {
                format!("Failed to write table to: {}", output_path.display())
            })?;

        writer.write_record(&self.headers)?;

        for row in &self.rows {
            writer.write_record(row)?;
        }

        Ok(())
    }
}

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
            .wrap_err(format!("Unable to write file: {:?}", tmp_path))?;
        decompress_file(&tmp_path, output_path, true)?;
    } else {
        let content = response.text().await?;
        write(output_path, content)
            .wrap_err(format!("Unable to write file: {:?}", output_path))?;
    }

    Ok(())
}

pub fn decompress_file(
    input: &PathBuf,
    output: &PathBuf,
    inplace: bool,
) -> Result<(), Report> {
    let ext = input.extension().unwrap();

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

/// Write table to file.
pub fn write_table(
    table: &Vec<Vec<String>>,
    output_path: &Path,
    delim: Option<char>,
) -> Result<(), Report> {
    let delim = match delim {
        Some(c) => c as u8,
        None => b'\t',
    };
    let mut writer = csv::WriterBuilder::new()
        .delimiter(delim)
        .from_path(output_path)
        .wrap_err_with(|| {
            format!("Failed to write table to: {}", output_path.display())
        })?;

    for row in table {
        writer.write_record(row)?;
    }

    Ok(())
}
