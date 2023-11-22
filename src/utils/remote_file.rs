use chrono::prelude::*;
use serde::{Deserialize, Serialize};
use std::default::Default;
use std::path::PathBuf;

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct RemoteFile {
    pub url: String,
    pub sha: String,
    pub local_path: PathBuf,
    pub date_created: DateTime<Utc>,
    pub date_downloaded: DateTime<Utc>,
}

impl Default for RemoteFile {
    fn default() -> Self {
        Self::new()
    }
}

impl RemoteFile {
    pub fn new() -> Self {
        RemoteFile {
            url: String::new(),
            sha: String::new(),
            local_path: PathBuf::new(),
            date_created: DateTime::default(),
            date_downloaded: DateTime::default(),
        }
    }
}
