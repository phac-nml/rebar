use crate::traits::ToYaml;
use color_eyre::eyre::Report;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::default::Default;

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum Mutation {
    Substitution,
    Deletion,
}

// ----------------------------------------------------------------------------
// Deletion
// ----------------------------------------------------------------------------

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Deletion {
    pub coord: usize,
    pub reference: char,
    pub alt: char,
}

impl std::fmt::Display for Deletion {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}{}{}", self.reference, self.coord, self.alt)
    }
}

impl PartialEq for Deletion {
    fn eq(&self, other: &Self) -> bool {
        self.coord == other.coord
            && self.reference == other.reference
            && self.alt == other.alt
    }
}

impl Eq for Deletion {}

impl Ord for Deletion {
    fn cmp(&self, other: &Self) -> Ordering {
        self.coord.cmp(&other.coord)
    }
}

impl PartialOrd for Deletion {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.coord.partial_cmp(&other.coord)
    }
}

// ----------------------------------------------------------------------------
// Substitution
// ----------------------------------------------------------------------------

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Substitution {
    pub coord: usize,
    pub reference: char,
    pub alt: char,
}

impl std::fmt::Display for Substitution {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}{}{}", self.reference, self.coord, self.alt)
    }
}

impl PartialEq for Substitution {
    fn eq(&self, other: &Self) -> bool {
        self.coord == other.coord
            && self.reference == other.reference
            && self.alt == other.alt
    }
}

impl Eq for Substitution {}

impl Ord for Substitution {
    fn cmp(&self, other: &Self) -> Ordering {
        self.coord.cmp(&other.coord)
    }
}

impl PartialOrd for Substitution {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.coord.partial_cmp(&other.coord)
    }
}

impl Substitution {
    pub fn to_deletion(&self) -> Deletion {
        Deletion {
            coord: self.coord,
            reference: self.reference,
            alt: '-',
        }
    }
}

// ----------------------------------------------------------------------------
// Substitution
// ----------------------------------------------------------------------------

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Sequence {
    pub id: String,
    pub seq: Vec<char>,
    alphabet: Vec<char>,
    pub substitutions: Vec<Substitution>,
    pub deletions: Vec<Deletion>,
    pub missing: Vec<usize>,
}

impl ToYaml for Sequence {
    fn to_yaml(&self) -> String {
        // This is terrible code formatting, I'm not surely how to improve yet,
        // While keeping str output format readable
        format!(
            "\n
            id:             {}
            missing:        {}
            deletions:      {}
            substitutions:  {}
            ",
            self.id,
            self.missing.iter().format(", "),
            self.deletions.iter().format(", "),
            self.substitutions.iter().format(", "),
        )
    }
}

impl Sequence {
    pub fn new() -> Self {
        Sequence {
            id: String::new(),
            seq: Vec::new(),
            alphabet: vec!['A', 'C', 'G', 'T'],
            substitutions: Vec::new(),
            deletions: Vec::new(),
            missing: Vec::new(),
        }
    }

    pub fn from_record(
        record: bio::io::fasta::Record,
        reference: Option<&Sequence>,
        mask: usize,
    ) -> Result<Self, Report> {
        let mut sample = Sequence::new();
        sample.id = record.id().to_string();
        sample.seq = record.seq().iter().map(|b| *b as char).collect();

        if let Some(reference) = reference {
            let genome_length = reference.seq.len();
            // Construct iterator to traverse sample and reference bases together
            let it = sample.seq.iter().zip(reference.seq.iter());
            for (i, (s, r)) in it.enumerate() {
                // Genomic coordinates are 1-based
                let coord: usize = i + 1;
                let mut s = *s;
                let r = *r;
                // Mask 5' and 3' ends
                if coord <= mask || coord > genome_length - mask {
                    s = 'N';
                }

                match s {
                    // Missing data (N)
                    'N' => sample.missing.push(coord),
                    // Reference Missing data (N)
                    _s if r == 'N' => continue,
                    // Deletion
                    '-' => {
                        let deletion = Deletion {
                            coord,
                            reference: r,
                            alt: s,
                        };
                        sample.deletions.push(deletion)
                    }
                    // Ambiguous data (IUPAC not in alphabet)
                    s if s != r && !sample.alphabet.contains(&s) => {
                        sample.missing.push(coord)
                    }
                    // Substitution
                    s if s != r => {
                        let substitution = Substitution {
                            coord,
                            reference: r,
                            alt: s,
                        };
                        sample.substitutions.push(substitution)
                    }
                    // Reference
                    _ => continue,
                }
            }
        }

        Ok(sample)
    }
}
