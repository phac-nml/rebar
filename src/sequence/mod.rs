pub mod parsimony;

use color_eyre::eyre::{eyre, ContextCompat, Report, Result, WrapErr};
use color_eyre::Help;
use noodles::{core::Position, fasta};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::default::Default;
use std::fs::File;
use std::io::BufReader;
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
        Some(self.cmp(other))
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
        Some(self.cmp(other))
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
        record: fasta::Record,
        reference: Option<&Sequence>,
        mask: &Vec<usize>,
    ) -> Result<Self, Report> {
        let mut sample = Sequence::new();

        // parse fasta::Record into rebar Sequence
        let id = record.name().to_string();

        // convert sequence to vec of bases, noodle positions are 1-based!
        let start = Position::try_from(1).unwrap();
        let seq: Vec<char> = record
            .sequence()
            .get(start..)
            .context(format!("Failed to parse sequence record {}", &sample.id))?
            .iter()
            .map(|b| *b as char)
            .collect();
        let sample_len = seq.len();
        sample.id = id.clone();
        sample.seq = seq;

        // check mask coord
        for bases in mask {
            if *bases > sample.seq.len() {
                return Err(
                    eyre!("5' and 3' masking ({mask:?}) is incompatible with {id} sequence length {sample_len}")
                    .suggestion("Please change your --mask parameter.")
                    .suggestion("Maybe you want to disable masking all together with --mask 0,0 ?")
                );
            }
        }

        if let Some(reference) = reference {
            let ref_len = reference.seq.len();
            if sample_len != ref_len {
                return Err(
                    eyre!("Reference and {id} are different lengths ({ref_len} vs {sample_len})!")
                    .suggestion(format!("Are you sure {id} is aligned correctly?"))
                );
            }
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
    let mut reader = File::open(path).map(BufReader::new).map(fasta::Reader::new)?;

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
