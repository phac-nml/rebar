use std::str::FromStr;

#[derive(Debug)]
pub enum LogVerbosity {
    Info,
    Warn,
    Debug,
}

impl std::fmt::Display for LogVerbosity {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Convert to lowercase for RUST_LOG env var compatibility
        let lowercase = format!("{:?}", self).to_lowercase();
        write!(f, "{lowercase}")
    }
}

#[derive(Debug)]
pub struct UnknownLogVerbosityError;

impl FromStr for LogVerbosity {
    type Err = UnknownLogVerbosityError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        match input {
            "debug" => Ok(LogVerbosity::Debug),
            "info" => Ok(LogVerbosity::Info),
            "warn" => Ok(LogVerbosity::Warn),
            _ => Err(UnknownLogVerbosityError),
        }
    }
}
