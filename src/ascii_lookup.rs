use std::num::ParseFloatError;

fn read_table<const NUM_ROWS: usize>(
    contents: &str,
) -> Result<Vec<[f64; NUM_ROWS]>, ParseFloatError> {
    let values: Result<Vec<_>, _> = contents.split_whitespace().map(|x| x.parse()).collect();
    let result = values?
        .chunks(NUM_ROWS)
        .map(|chunk| {
            let mut result = [0.0; NUM_ROWS];

            for i in 0..NUM_ROWS {
                result[i] = chunk[i]
            }
            result
        })
        .collect();

    Ok(result)
}

pub struct LookupTable<const NUM_ROWS: usize> {
    pub(crate) rows: Vec<[f64; NUM_ROWS]>,
}

impl<const NUM_ROWS: usize> LookupTable<NUM_ROWS> {
    /// Return a lookup table from a `Vec` of rows. `rows` must have length
    /// larger than 1, and the left-most column must be increasing
    /// monotonically.
    pub fn from_rows(rows: Vec<[f64; NUM_ROWS]>) -> Self {
        // ensure the table's left-column is ordered increasing monotonically.
        if rows.len() < 2 {
            panic!("the table must have at least two entries");
        }

        let mut x_prev = rows.first().unwrap()[0];
        for row in &rows[1..] {
            let x = row[0];
            if x <= x_prev {
                panic!("the table left column is not monotonically increasing");
            }
            x_prev = x;
        }
        Self { rows }
    }

    pub fn from_ascii(filename: &str) -> Self {
        let contents = std::fs::read_to_string(filename).unwrap();
        Self {
            rows: read_table(&contents).unwrap(),
        }
    }

    /// Return a fixed-length array of data at the given independent variable
    /// value `x`. The result is interpolated linearly between the two nearest
    /// tabulated points. This function panics if `x` is out of range (not
    /// between the lowest and highest value of the left-most table column.
    pub fn sample(&self, x: f64) -> [f64; NUM_ROWS] {
        let mut result = [0.0; NUM_ROWS];
        let (i0, i1) = self.indexes_straddling(x);
        let v = &self.rows;

        for i in 0..NUM_ROWS {
            let x0 = v[i0][0];
            let y0 = v[i0][i];
            let x1 = v[i1][0];
            let y1 = v[i1][i];
            result[i] = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
        }
        result
    }

    fn indexes_straddling(&self, x: f64) -> (usize, usize) {
        let xmin = self.rows.first().unwrap()[0];
        let xmax = self.rows.last().unwrap()[0];

        if x <= xmin {
            panic! {
                "attempt to sample table at or below smallest tabulated point ({} <= {})",
                x,
                xmin
            }
        }
        if x > xmax {
            panic! {
                "attempt to sample table above the largest tabulated point ({} > {})",
                x,
                xmax
            }
        }

        let index = match self
            .rows
            .binary_search_by(|row| Self::compare_f64(row[0], x))
        {
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
