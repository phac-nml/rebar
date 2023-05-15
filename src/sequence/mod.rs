use serde::{Serialize, Deserialize};
use crate::traits::ToYaml;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Sequence {
    pub id: String,
    pub seq: Vec<char>,  
}

impl ToYaml for Sequence {}