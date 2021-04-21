pub struct LookupTable {
    data: Vec<(f64, f64)>,
}

impl LookupTable {
    pub fn new(data: Vec<(f64, f64)>) -> Self {
        // ensure the table's left-column is ordered increasing monotonically.
        if data.len() < 2 {
            panic!("the table must have at least two enries");
        }

        let mut x_prev = data.first().unwrap().0;
        for &(x, _) in &data[1..] {
            if x <= x_prev {
                panic!("the table left column is not monotonically increasing");
            }
            x_prev = x;
        }
        Self { data }
    }

    pub fn sample(&self, x: f64) -> f64 {
        let (i0, i1) = self.indexes_straddling(x);
        let v = &self.data;
        let x0 = v[i0].0;
        let y0 = v[i0].1;
        let x1 = v[i1].0;
        let y1 = v[i1].1;
        y0 + (x - x0) * (y1 - y0) / (x1 - x0)
    }

    fn indexes_straddling(&self, x: f64) -> (usize, usize) {
        if x <= self.data[0].0 {
            panic! {
                "attempt to sample table at or below smallest tabulated point ({} <= {})",
                x,
                self.data[0].0
            };
        }

        let index = match self.data.binary_search_by(|&(xi, _)| Self::compare_f64(xi, x)) {
            Ok(index) => index,
            Err(index) => index,
        };

        (index - 1, index)
    }

    fn compare_f64(a: f64, b: f64) -> std::cmp::Ordering {
        if a < b {
            std::cmp::Ordering::Less
        } else if a > b {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn lookup_table_panics_if_sampled_at_lower_bound() {
        let table = LookupTable::new(vec![(0.0, 0.1), (1.0, 0.2), (2.0, 0.3)]);
        table.indexes_straddling(0.0);
    }

    #[test]
    fn lookup_table_does_not_panic_if_sampled_at_upper_bound() {
        let table = LookupTable::new(vec![(0.0, 0.1), (1.0, 0.2), (2.0, 0.3)]);
        table.indexes_straddling(2.0);
    }

    #[test]
    fn lookup_table_gives_the_right_indexes_straddling() {
        let table = LookupTable::new(vec![(0.0, 0.1), (1.0, 0.2), (2.0, 0.3)]);
        assert_eq!(table.indexes_straddling(0.5), (0, 1));
        assert_eq!(table.indexes_straddling(1.0), (0, 1));
        assert_eq!(table.indexes_straddling(1.5), (1, 2));
    }

    #[test]
    fn lookup_table_can_be_sampled_at_tabulated_points() {
        let table = LookupTable::new(vec![(0.0, 0.1), (1.0, 0.2), (2.0, 0.3)]);
        assert!(f64::abs(table.sample(0.5) - 0.15) < 1e-10);
        assert!(f64::abs(table.sample(1.0) - 0.20) < 1e-10);
        assert!(f64::abs(table.sample(1.5) - 0.25) < 1e-10);
    }
}
