use std::path::Path;
use serde::{Serialize, Deserialize};




// ============================================================================
#[derive(thiserror::Error, Debug)]
pub enum Error {

    #[error("{0}")]
    CiboriumSer(#[from] ciborium::ser::Error<std::io::Error>),

    #[error("{0}")]
    CiboriumDe(#[from] ciborium::de::Error<std::io::Error>),

    #[error("{0}")]
    IO(#[from] std::io::Error),
}




// ============================================================================
pub fn parent_directory(path_str: &str) -> String {
    match Path::new(&path_str).parent().and_then(Path::to_str) {
        None     => ".",
        Some("") => ".",
        Some(parent) => parent,
    }.into()
}

pub fn write_cbor<T: Serialize>(value: &T, path_str: &str) -> Result<(), Error> {
    println!("write {}", path_str);
    let file = std::fs::File::create(&path_str)?;
    let mut buffer = std::io::BufWriter::new(file);
    Ok(ciborium::ser::into_writer(&value, &mut buffer)?)
}

pub fn read_cbor<T: for<'de> Deserialize<'de>>(path_str: &str) -> Result<T, Error> {
    let file = std::fs::File::open(path_str)?;
    let buffer = std::io::BufReader::new(file);
    Ok(ciborium::de::from_reader(buffer)?)
}
