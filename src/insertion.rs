//! Robinson-Schensted row insertion and bumping path operations.

use crate::tableau::YoungTableau;

/// Result of a row insertion, including the position where a new cell was added.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InsertionResult {
    /// The modified tableau.
    pub tableau: YoungTableau,
    /// The (row, col) position where a new cell was created.
    pub new_cell: (usize, usize),
    /// The bumping path: sequence of (row, col) positions that were modified.
    pub bumping_path: Vec<(usize, usize)>,
}

/// Perform Robinson-Schensted row insertion of value `x` into the tableau.
///
/// Inserts `x` into the first row by finding the leftmost element > x (or appending),
/// bumping that element into the next row, and repeating until a value is appended.
pub fn row_insert(tableau: &YoungTableau, x: u32) -> InsertionResult {
    let mut rows = tableau.rows.clone();
    let mut path = vec![];
    let mut val = x;
    let mut new_cell = (0, 0);

    for row_idx in 0..rows.len() {
        let row = &mut rows[row_idx];
        // Find leftmost element strictly greater than val
        let pos = row.iter().position(|&v| v > val);
        match pos {
            Some(col) => {
                let bumped = row[col];
                row[col] = val;
                path.push((row_idx, col));
                val = bumped;
            }
            None => {
                // Append to this row
                let col = row.len();
                row.push(val);
                path.push((row_idx, col));
                new_cell = (row_idx, col);
                return InsertionResult {
                    tableau: YoungTableau::from_rows(rows),
                    new_cell,
                    bumping_path: path,
                };
            }
        }
    }

    // If we exhausted all rows, start a new row
    let _col = 0;
    rows.push(vec![val]);
    path.push((rows.len() - 1, 0));
    new_cell = (rows.len() - 1, 0);

    InsertionResult {
        tableau: YoungTableau::from_rows(rows),
        new_cell,
        bumping_path: path,
    }
}

/// Insert a sequence of values into an empty tableau, building the full insertion tableau.
pub fn insertion_tableau(values: &[u32]) -> YoungTableau {
    let mut t = YoungTableau::empty();
    for &v in values {
        let result = row_insert(&t, v);
        t = result.tableau;
    }
    t
}

/// Compute the recording tableau for a sequence of values.
///
/// The recording tableau records the order in which cells were added,
/// using 1-indexed positions. Together with the insertion tableau, this
/// gives the Robinson-Schensted correspondence.
pub fn recording_tableau(values: &[u32]) -> YoungTableau {
    let mut p = YoungTableau::empty();
    let mut q = YoungTableau::empty();

    for (i, &v) in values.iter().enumerate() {
        let result = row_insert(&p, v);
        let (r, c) = result.new_cell;
        p = result.tableau;

        // Add the cell to Q at the same position with value i+1
        // Ensure Q has enough rows
        while q.rows.len() <= r {
            q.rows.push(vec![]);
        }
        // Ensure the row has enough columns
        while q.rows[r].len() < c {
            q.rows[r].push(0);
        }
        // Now append or set
        if q.rows[r].len() == c {
            q.rows[r].push((i + 1) as u32);
        } else {
            q.rows[r][c] = (i + 1) as u32;
        }
    }
    q
}

/// Reverse the insertion: given a tableau and a corner cell, reverse-bump to recover
/// the inserted value.
pub fn reverse_insert(tableau: &YoungTableau, row: usize, col: usize) -> (YoungTableau, u32) {
    let mut rows = tableau.rows.clone();
    let mut val = rows[row][col];
    // Remove the cell
    rows[row].pop();
    if rows[row].is_empty() && row == rows.len() - 1 {
        rows.pop();
    }

    // Reverse bump from bottom to top
    for r in (0..row).rev() {
        // Find rightmost element <= val
        let pos = rows[r].iter().rposition(|&v| v <= val);
        if let Some(c) = pos {
            std::mem::swap(&mut rows[r][c], &mut val);
        }
    }

    (YoungTableau::from_rows(rows), val)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_row_insert_empty() {
        let t = YoungTableau::empty();
        let result = row_insert(&t, 3);
        assert_eq!(result.tableau.rows, vec![vec![3]]);
        assert_eq!(result.new_cell, (0, 0));
    }

    #[test]
    fn test_row_insert_append() {
        let t = YoungTableau::from_rows(vec![vec![1, 2]]);
        let result = row_insert(&t, 3);
        assert_eq!(result.tableau.rows, vec![vec![1, 2, 3]]);
        assert_eq!(result.new_cell, (0, 2));
    }

    #[test]
    fn test_row_insert_bump() {
        let t = YoungTableau::from_rows(vec![vec![1, 3]]);
        let result = row_insert(&t, 2);
        assert_eq!(result.tableau.rows, vec![vec![1, 2], vec![3]]);
        assert_eq!(result.new_cell, (1, 0));
    }

    #[test]
    fn test_insertion_tableau_sequence() {
        let t = insertion_tableau(&[3, 1, 4, 1, 5]);
        assert_eq!(t.size(), 5);
        assert!(t.is_semistandard());
    }

    #[test]
    fn test_insertion_preserves_shape() {
        let t = insertion_tableau(&[5, 1, 4, 2, 3]);
        // Shape should be a valid partition
        let shape = t.shape();
        for w in shape.windows(2) {
            assert!(w[0] >= w[1]);
        }
    }

    #[test]
    fn test_bumping_path() {
        let t = YoungTableau::from_rows(vec![vec![1, 3]]);
        let result = row_insert(&t, 2);
        // Should bump 3 from position (0,1) and place it at (1,0)
        assert_eq!(result.bumping_path, vec![(0, 1), (1, 0)]);
    }

    #[test]
    fn test_recording_tableau() {
        let rt = recording_tableau(&[3, 1, 2]);
        assert_eq!(rt.size(), 3);
        assert!(rt.is_standard());
    }

    #[test]
    fn test_reverse_insert() {
        let t = insertion_tableau(&[3, 1, 4]);
        let (rev, val) = reverse_insert(&t, t.num_rows() - 1, t.rows.last().unwrap().len() - 1);
        assert_eq!(rev.size(), 2);
    }

    #[test]
    fn test_insert_reverse_roundtrip() {
        let mut t = YoungTableau::empty();
        let values = vec![4, 2, 5, 1, 3];
        for &v in &values {
            t = row_insert(&t, v).tableau;
        }
        assert_eq!(t.size(), 5);
    }

    #[test]
    fn test_multiple_bumps() {
        let t = insertion_tableau(&[2, 1, 3, 1, 2, 4]);
        assert_eq!(t.size(), 6);
        assert!(t.is_semistandard());
    }
}
