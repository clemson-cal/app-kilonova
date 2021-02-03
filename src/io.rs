use std::path::Path;
use serde::{Serialize, Deserialize};




// ============================================================================
#[derive(thiserror::Error, Debug)]
pub enum Error {

    #[error("{0}")]
    SerdeCbor(#[from] serde_cbor::Error),

    #[error("{0}")]
    IO(#[from] std::io::Error),

    #[error("input file is compressed, but snap is not enabled")]
    CannotReadSnappy,

    #[error("snappy_compression = true, but snap is not enabled")]
    CannotWriteSnappy,
}




// ============================================================================
pub fn parent_directory(path_str: &str) -> String {
    match Path::new(&path_str).parent().and_then(Path::to_str) {
        None     => ".",
        Some("") => ".",
        Some(parent) => parent,
    }.into()
}

#[cfg(feature = "serde_cbor")]
pub fn write_cbor<T: Serialize>(value: &T, path_str: &str, snappy_compression: bool) -> Result<(), Error> {
    println!("write {}", path_str);
    let file = std::fs::File::create(&path_str)?;
    let buffer = std::io::BufWriter::new(file);

    if snappy_compression {
        #[cfg(feature = "snap")] {
            serde_cbor::to_writer(snap::write::FrameEncoder::new(buffer), &value)?;
        }
        #[cfg(not(feature = "snap"))] {
            return Err(Error::CannotWriteSnappy)
        }
    } else {
        serde_cbor::to_writer(buffer, &value)?;        
    }
    Ok(())
}

#[cfg(not(feature = "serde_cbor"))]
pub fn write_cbor<T: Serialize>(_: &T, path_str: &str, _: bool) -> Result<(), Error> {
    println!("skip writing {} (serde_cbor is not enabled)", path_str);
    Ok(())
}

#[cfg(feature = "serde_cbor")]
pub fn read_cbor<T: for<'de> Deserialize<'de>>(path_str: &str, snappy_compression: bool) -> Result<T, Error> {
    let file = std::fs::File::open(path_str)?;
    let buffer = std::io::BufReader::new(file);

    if snappy_compression {
        #[cfg(feature = "snap")] {
            Ok(serde_cbor::from_reader(snap::read::FrameDecoder::new(buffer))?)
        }
        #[cfg(not(feature = "snap"))] {
            return Err(Error::CannotReadSnappy)
        }
    } else {
        Ok(serde_cbor::from_reader(buffer)?)
    }
}

#[cfg(not(feature = "serde_cbor"))]
pub fn read_cbor<T: for<'de> Deserialize<'de>>(path_str: &str, _: bool) -> Result<T, Error> {
    anyhow::bail!("input file {} given, but serde_cbor is not enabled", path_str)
}
