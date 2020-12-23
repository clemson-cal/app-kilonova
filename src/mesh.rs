use std::cmp::PartialEq;
use std::str::FromStr;
use std::fmt::{Display, Formatter};




// ============================================================================
#[derive(Copy, Clone, Eq, Hash)]
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

impl Display for BlockIndex {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "0:{:03}-{:03}", self.i, self.j)
    }
}

impl FromStr for BlockIndex {
    type Err = <usize as FromStr>::Err;

    fn from_str(string: &str) -> std::result::Result<Self, Self::Err> {
        // 0:000-000
        let i = usize::from_str(&string[2..5])?;
        let j = usize::from_str(&string[6..9])?;
        Ok(Self{i, j})
    }
}




// ============================================================================
mod test {
    use crate::mesh::BlockIndex;
    use std::str::FromStr;

    #[test]
    fn block_index_can_be_parsed() {
        assert!(BlockIndex::from_str("0:000-000").unwrap() == BlockIndex{i:0, j:0});
        assert!(BlockIndex::from_str("0:666-012").unwrap() == BlockIndex{i:666, j:12});
    }

    #[test]
    fn block_index_goes_to_string() {
        assert!(BlockIndex{i:666, j:12}.to_string() == "0:666-012");
    }

    #[test]
    #[should_panic]
    fn bad_block_index_panics() {
        BlockIndex::from_str("0:000-").unwrap();
    }
}
