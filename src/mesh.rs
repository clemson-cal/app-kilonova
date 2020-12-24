use std::cmp::PartialEq;
use std::str::FromStr;
use std::fmt::{Display, Formatter};
use serde::{Serialize, Deserialize};




// ============================================================================
#[derive(Copy, Clone, Eq, Hash, Serialize, Deserialize)]
pub struct BlockIndex {
    i: usize,
    j: usize,
}




// ============================================================================
impl PartialEq for BlockIndex {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i && self.j == other.j
    }
}
