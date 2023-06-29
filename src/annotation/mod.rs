use crate::dataset::Name;
use crate::utils;
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;

#[derive(Clone, Debug)]
pub struct Annotations {
    gene: Vec<String>,
    abbreviation: Vec<String>,
    start: Vec<usize>,
    end: Vec<usize>,
}

impl Annotations {
    pub fn from_name(name: &Name) -> Result<Self, Report> {
        let annotations = match name {
            // ----------------------------------------------------------------
            // sars-cov-2 annotations
            Name::SarsCov2 => Annotations {
                gene: vec![
                    "ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b",
                    "ORF8", "ORF9b",
                ]
                .into_iter()
                .map(String::from)
                .collect_vec(),
                abbreviation: vec![
                    "1a", "1b", "S", "3a", "E", "M", "6", "7a", "7b", "8", "9b",
                ]
                .into_iter()
                .map(String::from)
                .collect_vec(),
                start: vec![
                    266, 13468, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894,
                    28284,
                ],
                end: vec![
                    13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259,
                    28577,
                ],
            },
            _ => {
                return Err(eyre!(
                    "Annotations for the {name} dataset is not implemented yet."
                ))
            }
        };

        Ok(annotations)
    }

    pub fn to_table(&self) -> Result<utils::Table, Report> {
        let mut table = utils::Table::new();

        // headers
        table.headers = vec!["gene", "abbreviation", "start", "end"]
            .into_iter()
            .map(String::from)
            .collect_vec();

        // rows
        let mut rows = Vec::new();
        for i in 0..(self.gene.len()) {
            let row = vec![
                self.gene[i].clone(),
                self.abbreviation[i].clone(),
                self.start[i].to_string(),
                self.end[i].to_string(),
            ];
            rows.push(row);
        }

        table.rows = rows;

        Ok(table)
    }
}
