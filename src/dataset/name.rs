use color_eyre::eyre::{Report, Result};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
pub enum Name {
    SarsCov2,
    RsvA,
    RsvB,
    Unknown,
}

impl fmt::Display for Name {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Name::SarsCov2 => String::from("sars-cov-2"),
            Name::RsvA => String::from("rsv-a"),
            Name::RsvB => String::from("rsv-b"),
            Name::Unknown => String::from("unknown"),
        };

        write!(f, "{}", name)
    }
}

impl FromStr for Name {
    type Err = Report;

    fn from_str(name: &str) -> Result<Self, Report> {
        let name = match name {
            "sars-cov-2" => Name::SarsCov2,
            "rsv-a" => Name::RsvA,
            "rsv-b" => Name::RsvB,
            _ => Name::Unknown,
        };

        Ok(name)
    }
}
