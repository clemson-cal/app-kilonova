use std::io::Read;
use serde::{Serialize, Deserialize};
use serde_yaml::{Value, Mapping, to_value, from_value, from_str, from_reader};




// ============================================================================
fn merge_mapping(value_map: &Mapping, patch_map: &Mapping) -> Mapping {
    let mut result = value_map.clone();

    for (key, patch_value) in patch_map {
        let new_value = merge_value(value_map.get(key).unwrap_or(&Value::Null), patch_value);
        result.insert(key.clone(), new_value);
    }
    return result;
}




// ============================================================================
fn merge_value(value: &Value, patch: &Value) -> Value {
    if let (Some(value_map), Some(patch_map)) = (value.as_mapping(), patch.as_mapping()) {
        Value::from(merge_mapping(value_map, patch_map))
    } else {
        patch.clone()
    }
}




/**
 * Extends anything that is Clone, Serialize, and Deserialize to have mutable
 * "patch" methods, accepting `serde_yaml::Value` objects.
 */
pub trait Patch {
    fn patch_from_value(&mut self, patch_value: &Value) -> serde_yaml::Result<()>;
    fn patch_from_str(&mut self, yaml_str: &str) -> serde_yaml::Result<()> {
        self.patch_from_value(&from_str(yaml_str)?)
    }
    fn patch_from_reader<R>(&mut self, reader: R) -> serde_yaml::Result<()> where R: Read {
        self.patch_from_value(&from_reader(reader)?)
    }
}




// ============================================================================
impl<T> Patch for T where T: Clone + Serialize + for<'de> Deserialize<'de> {
    fn patch_from_value(&mut self, patch_value: &Value) -> serde_yaml::Result<()> {
        let self_value = to_value(self.clone())?;
        let merged_self_value = merge_value(&self_value, &patch_value);
        let merged_self: T = from_value(merged_self_value)?;
        *self = merged_self;
        Ok(())
    }
}




// ============================================================================
mod test {
    #[derive(Clone, serde::Serialize, serde::Deserialize)]
    struct Config {
        x: f64,
        y: usize,
    }

    #[test]
    fn can_merge() {
        use crate::yaml_patch::Patch;
        let mut config = Config{x: 32.0, y: 512};
        config.patch_from_str(r"y: 1024").unwrap();
        assert!(config.x == 32.0);
        assert!(config.y == 1024);
    }
}
