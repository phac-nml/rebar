use color_eyre::eyre::{eyre, Report, Result};
use std::str::FromStr;

#[derive(Debug)]
pub enum Verbosity {
    Info,
    Warn,
    Debug,
}

impl std::fmt::Display for Verbosity {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Convert to lowercase for RUST_LOG env var compatibility
        let lowercase = format!("{:?}", self).to_lowercase();
        write!(f, "{lowercase}")
    }
}

impl FromStr for Verbosity {
    type Err = Report;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        match input {
            "debug" => Ok(Verbosity::Debug),
            "info" => Ok(Verbosity::Info),
            "warn" => Ok(Verbosity::Warn),
            _ => Err(eyre!("Unknown verbosity level {input}.")),
        }
    }
}
