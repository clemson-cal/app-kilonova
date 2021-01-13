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

fn merge_value(value: &Value, patch: &Value) -> Value {
    if let (Some(value_map), Some(patch_map)) = (value.as_mapping(), patch.as_mapping()) {
        Value::from(merge_mapping(value_map, patch_map))
    } else {
        patch.clone()
    }
}

fn single_item_mapping<T: Into<Value>>(val: T, key: &str) -> Mapping {
    let mut result = Mapping::new();
    result.insert(Value::from(key), val.into());
    result
}




/**
 * Extends a type with mutable key-value symantics, allowing it to be "patched"
 * with runtime yaml-derived values.
 */
pub trait Patch {

    /**
     * Update this object to reflect the contents of a `serde_yaml::Value`. The
     * rules are as follows:
     *
     * - If `patch_value` is not a `serde_yaml::Mapping`, then it is
     * deserialized into self.
     *
     * - Otherwise, if `patch_value` is a `serde_yaml::Mapping`, then each of
     * its values are merged into their respective data members of `self` using
     * this algorithm recursively.
     */
    fn patch_from_value(&mut self, patch_value: &Value) -> serde_yaml::Result<()>;

    /**
     * Update this object to reflect the contents of a YAML-encoded string.
     */
    fn patch_from_str(&mut self, yaml_str: &str) -> serde_yaml::Result<()> {
        self.patch_from_value(&from_str(yaml_str)?)
    }

    /**
     * Update this object to reflect the contents of a reader, which will yield
     * a YAML-encoded string.
     */
    fn patch_from_reader<R>(&mut self, reader: R) -> serde_yaml::Result<()> where R: Read {
        self.patch_from_value(&from_reader(reader)?)
    }

    /**
     * Update a single, possibly nested, data member within this object using
     * key-path style attribute access. For example, the string
     * "company.ceo.name = Bob" would operate on a data member called `company`
     * within this object, setting the `name` field of the `ceo` field to
     * `"Bob"`.
     */
    fn patch_from_key_val(&mut self, key_val_str: &str) -> anyhow::Result<()> {
        let tokens: Vec<_> = key_val_str.split('=').collect();

        if tokens.len() != 2 || tokens[1].is_empty() {
            anyhow::bail!("badly formed key=value in '{}'", key_val_str)
        }

        let mut parts = tokens[0].split('.').rev();
        let value: Value = from_str(tokens[1])?;
        let mapping = single_item_mapping(value, parts.next().unwrap());

        Ok(self.patch_from_value(&parts.fold(mapping, single_item_mapping).into())?)
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
    struct Inner {
        a: String,
        b: String,
    }

    #[derive(Clone, serde::Serialize, serde::Deserialize)]
    struct Config {
        x: f64,
        y: usize,
        inner: Inner,
    }

    impl Default for Config {
        fn default() -> Self {
            Config{x: 32.0, y: 512, inner: Inner{a: "a".into(), b: "b".into()}}
        }
    }

    #[test]
    fn can_merge_from_str() {
        use crate::yaml_patch::Patch;
        let mut config = Config::default();
        config.patch_from_str(r"y: 1024").unwrap();
        assert!(config.x == 32.0);
        assert!(config.y == 1024);
        assert!(config.inner.a == "a");
        assert!(config.inner.b == "b");
    }

    #[test]
    fn can_merge_from_key_val() {
        use crate::yaml_patch::Patch;
        let mut config = Config::default();
        config.patch_from_key_val(r"y=1024").unwrap();
        config.patch_from_key_val(r"inner.a=A").unwrap();
        assert!(config.x == 32.0);
        assert!(config.y == 1024);
    }
}
