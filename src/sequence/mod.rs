use crate::traits::ToYaml;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Sequence {
    pub id: String,
    pub seq: Vec<char>,
}

impl ToYaml for Sequence {}
