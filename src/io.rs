use std::path::Path;
use serde::{Serialize, Deserialize};




// ============================================================================
pub fn parent_directory(path_str: &str) -> String {
    match Path::new(&path_str).parent().and_then(Path::to_str) {
        None     => ".",
        Some("") => ".",
        Some(parent) => parent,
    }.into()
}

#[cfg(feature = "serde_cbor")]
pub fn write_cbor<T: Serialize>(value: &T, path_str: &str, snappy_compression: bool) -> anyhow::Result<()> {
    println!("write {}", path_str);
    let file = std::fs::File::create(&path_str)?;
    let buffer = std::io::BufWriter::new(file);

    if snappy_compression {
        #[cfg(feature = "snap")] {
            serde_cbor::to_writer(snap::write::FrameEncoder::new(buffer), &value)?;
        }
        #[cfg(not(feature = "snap"))] {
            anyhow::bail!("snappy_compression = true, but snap is not enabled")
        }
    } else {
        serde_cbor::to_writer(buffer, &value)?;        
    }
    Ok(())
}

#[cfg(not(feature = "serde_cbor"))]
pub fn write_cbor<T: Serialize>(_: &T, path_str: &str, _: bool) -> anyhow::Result<()> {
    println!("skip writing {} (serde_cbor is not enabled)", path_str);
    Ok(())
}

#[cfg(feature = "serde_cbor")]
pub fn read_cbor<T: for<'de> Deserialize<'de>>(path_str: &str, snappy_compression: bool) -> anyhow::Result<T> {
    let file = std::fs::File::open(path_str)?;
    let buffer = std::io::BufReader::new(file);

    if snappy_compression {
        #[cfg(feature = "snap")] {
            Ok(serde_cbor::from_reader(snap::read::FrameDecoder::new(buffer))?)
        }
        #[cfg(not(feature = "snap"))] {
            anyhow::bail!("input file is compressed, but snap is not enabled")
        }
    } else {
        Ok(serde_cbor::from_reader(buffer)?)
    }
}

#[cfg(not(feature = "serde_cbor"))]
pub fn read_cbor<T: for<'de> Deserialize<'de>>(path_str: &str, _: bool) -> anyhow::Result<T> {
    anyhow::bail!("input file {} given, but serde_cbor is not enabled", path_str)
}
