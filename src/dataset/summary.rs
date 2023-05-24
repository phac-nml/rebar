use crate::traits::ToYaml;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Summary {
    pub tag: String,
    pub name: String,
}

impl ToYaml for Summary {}
