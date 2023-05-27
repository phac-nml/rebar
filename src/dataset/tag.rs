use color_eyre::eyre::{Report, Result};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub enum Tag {
    Latest,
    Nightly,
    Archive(String),
    Unknown,
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let tag = match self {
            Tag::Latest => String::from("latest"),
            Tag::Nightly => String::from("nightly"),
            Tag::Archive(tag) => tag.to_owned(),
            Tag::Unknown => String::from("unknown"),
        };

        write!(f, "{}", tag)
    }
}

impl FromStr for Tag {
    type Err = Report;

    fn from_str(tag: &str) -> Result<Tag, Report> {
        let tag = match tag {
            "latest" => Tag::Latest,
            "nightly" => Tag::Nightly,
            _ => Tag::Archive(String::from(tag)),
        };

        Ok(tag)
    }
}
