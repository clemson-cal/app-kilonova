use crate::traits;

impl traits::Zeros for hydro_srhd::srhd_2d::Conserved {
    fn zeros() -> Self {
    	Self::default()
    }
}

impl traits::Arithmetic for hydro_srhd::srhd_2d::Conserved {
}

impl traits::Conserved for hydro_srhd::srhd_2d::Conserved {
    fn mass_and_momentum(&self) -> (f64, f64, f64) {
    	todo!()
    }
}
