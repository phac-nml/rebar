use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::default::Default;
use std::path::Path;
use std::str::FromStr;

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

#[derive(Copy, Clone, Debug, Hash, Serialize, Deserialize, PartialEq)]
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

impl FromStr for Substitution {
    type Err = Report;

    fn from_str(text: &str) -> Result<Self, Report> {
        let reference = text.chars().next().unwrap();
        let alt = text.chars().nth(text.len() - 1).unwrap();
        let coord = text[1..text.len() - 1].parse().unwrap();
        let substitution = Substitution {
            reference,
            alt,
            coord,
        };

        Ok(substitution)
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
    pub genome_length: usize,
    pub substitutions: Vec<Substitution>,
    pub deletions: Vec<Deletion>,
    pub missing: Vec<usize>,
}

impl Sequence {
    pub fn new() -> Self {
        Sequence {
            id: String::new(),
            seq: Vec::new(),
            alphabet: vec!['A', 'C', 'G', 'T'],
            genome_length: 0,
            substitutions: Vec::new(),
            deletions: Vec::new(),
            missing: Vec::new(),
        }
    }

    pub fn from_record(
        record: bio::io::fasta::Record,
        reference: Option<&Sequence>,
        mask: &Vec<usize>,
    ) -> Result<Self, Report> {
        let mut sample = Sequence::new();
        sample.id = record.id().to_string();
        sample.seq = record.seq().iter().map(|b| *b as char).collect();

        if let Some(reference) = reference {
            sample.genome_length = reference.seq.len();
            // Construct iterator to traverse sample and reference bases together
            let it = sample.seq.iter().zip(reference.seq.iter());
            for (i, (s, r)) in it.enumerate() {
                // Genomic coordinates are 1-based
                let coord: usize = i + 1;
                let mut s = *s;
                let r = *r;
                // Mask 5' and 3' ends
                if !mask.is_empty() && coord <= mask[0] {
                    s = 'N';
                }
                if mask.len() == 2 && coord > sample.genome_length - mask[1] {
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
        } else {
            sample.genome_length = sample.seq.len();
        }

        Ok(sample)
    }
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

/// Read first record of fasta path into sequence record.
pub fn read_reference(path: &Path, mask: &Vec<usize>) -> Result<Sequence, Report> {
    // start reading in the reference as fasta, raise error if file doesn't exist
    let reader = fasta::Reader::from_file(path).expect("Unable to read reference");

    // parse just the first record from the reference
    // 1. raise error if record iterator doesn't work
    // 2. raise error if first record is not proper fasta format.
    let reference = reader
        .records()
        .next()
        .ok_or_else(|| eyre!("Unable to read reference records: {path:?}"))?
        .wrap_err_with(|| eyre!("Unable to read first fasta record: {path:?}"))?;

    // convert to sequence
    let reference = Sequence::from_record(reference, None, mask)?;

    Ok(reference)
}
