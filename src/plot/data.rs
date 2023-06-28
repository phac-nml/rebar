use color_eyre::eyre::{eyre, Report, Result};
use csv;
use itertools::Itertools;
use std::path::{Path, PathBuf};

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
        let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;

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
}
