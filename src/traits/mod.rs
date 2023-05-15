use serde::Serialize;

pub trait ToYaml {
    fn to_yaml(&self) -> String
    where
        Self: Serialize,
    {
        serde_yaml::to_string(&self).unwrap()
    }
}
