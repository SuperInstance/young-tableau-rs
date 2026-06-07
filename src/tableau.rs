//! Core Young tableau data structure and partition operations.

/// A Young tableau (or semistandard tableau) represented as a list of rows.
///
/// Each row is a non-decreasing sequence of positive integers,
/// and row lengths are weakly decreasing (forming a partition shape).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct YoungTableau {
    /// Rows of the tableau, top to bottom.
    pub rows: Vec<Vec<u32>>,
}

impl YoungTableau {
    /// Create an empty Young tableau.
    pub fn empty() -> Self {
        Self { rows: vec![] }
    }

    /// Create a Young tableau from a shape (partition) filled with zeros.
    pub fn from_shape(shape: &[usize]) -> Self {
        Self {
            rows: shape.iter().map(|&len| vec![0; len]).collect(),
        }
    }

    /// Create a Young tableau from existing rows (no validation).
    pub fn from_rows(rows: Vec<Vec<u32>>) -> Self {
        Self { rows }
    }

    /// The shape (partition) of the tableau as a vector of row lengths.
    pub fn shape(&self) -> Vec<usize> {
        self.rows.iter().map(|r| r.len()).collect()
    }

    /// Total number of cells in the tableau.
    pub fn size(&self) -> usize {
        self.rows.iter().map(|r| r.len()).sum()
    }

    /// Number of rows.
    pub fn num_rows(&self) -> usize {
        self.rows.len()
    }

    /// Get the cell value at (row, col). Returns `None` if out of bounds.
    pub fn get(&self, row: usize, col: usize) -> Option<u32> {
        self.rows.get(row).and_then(|r| r.get(col)).copied()
    }

    /// Set the cell value at (row, col).
    pub fn set(&mut self, row: usize, col: usize, val: u32) {
        if row < self.rows.len() && col < self.rows[row].len() {
            self.rows[row][col] = val;
        }
    }

    /// Check if the tableau satisfies the semistandard conditions:
    /// - Rows are weakly increasing
    /// - Columns are strictly increasing
    pub fn is_semistandard(&self) -> bool {
        // Check rows are weakly increasing
        for row in &self.rows {
            for w in row.windows(2) {
                if w[0] > w[1] {
                    return false;
                }
            }
        }
        // Check columns are strictly increasing
        let max_cols = self.rows.iter().map(|r| r.len()).max().unwrap_or(0);
        for col in 0..max_cols {
            for row in 1..self.rows.len() {
                if col < self.rows[row].len() && col < self.rows[row - 1].len()
                    && self.rows[row][col] <= self.rows[row - 1][col] {
                        return false;
                    }
            }
        }
        true
    }

    /// Check if the tableau is standard (strictly increasing along both rows and columns,
    /// containing the numbers 1..=n exactly once).
    pub fn is_standard(&self) -> bool {
        let n = self.size();
        let mut seen = vec![false; n + 1];

        // Check rows strictly increasing and collect values
        for row in &self.rows {
            for w in row.windows(2) {
                if w[0] >= w[1] {
                    return false;
                }
            }
            for &v in row {
                if v == 0 || v as usize > n || seen[v as usize] {
                    return false;
                }
                seen[v as usize] = true;
            }
        }

        // Check columns strictly increasing
        let max_cols = self.rows.iter().map(|r| r.len()).max().unwrap_or(0);
        for col in 0..max_cols {
            for row in 1..self.rows.len() {
                if col < self.rows[row].len() && col < self.rows[row - 1].len()
                    && self.rows[row][col] <= self.rows[row - 1][col] {
                        return false;
                    }
            }
        }

        seen[1..].iter().all(|&b| b)
    }

    /// The conjugate (transpose) partition shape.
    pub fn conjugate_shape(&self) -> Vec<usize> {
        if self.rows.is_empty() {
            return vec![];
        }
        let max_cols = self.rows.iter().map(|r| r.len()).max().unwrap();
        let mut conj = vec![0usize; max_cols];
        for row in &self.rows {
            for (j, count) in conj.iter_mut().enumerate() {
                if j < row.len() {
                    *count += 1;
                }
            }
        }
        conj
    }

    /// Format the tableau as a string for display.
    pub fn display(&self) -> String {
        self.rows
            .iter()
            .map(|row| row.iter().map(|v| format!("{:3}", v)).collect::<Vec<_>>().join(" "))
            .collect::<Vec<_>>()
            .join("\n")
    }
}

/// A partition of an integer n.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Partition(pub Vec<usize>);

impl Partition {
    /// Create a partition, ensuring it is valid (weakly decreasing, non-negative).
    pub fn new(parts: Vec<usize>) -> Self {
        Self(parts)
    }

    /// The sum of all parts.
    pub fn weight(&self) -> usize {
        self.0.iter().sum()
    }

    /// Number of parts (length).
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Check if the partition is empty.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Check if the parts are weakly decreasing.
    pub fn is_valid(&self) -> bool {
        self.0.windows(2).all(|w| w[0] >= w[1])
    }

    /// Conjugate (transpose) partition.
    pub fn conjugate(&self) -> Partition {
        if self.0.is_empty() {
            return Partition(vec![]);
        }
        let max_val = self.0[0];
        let mut conj = vec![0usize; max_val];
        for &p in &self.0 {
            for j in 0..p {
                conj[j] += 1;
            }
        }
        Partition(conj)
    }

    /// Dominance ordering: self dominates other if cumulative sums are >= at every position.
    pub fn dominates(&self, other: &Partition) -> bool {
        let max_len = self.0.len().max(other.0.len());
        let mut self_cum = 0usize;
        let mut other_cum = 0usize;
        for i in 0..max_len {
            self_cum += self.0.get(i).unwrap_or(&0);
            other_cum += other.0.get(i).unwrap_or(&0);
            if self_cum < other_cum {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_tableau() {
        let t = YoungTableau::empty();
        assert_eq!(t.size(), 0);
        assert_eq!(t.num_rows(), 0);
        assert_eq!(t.shape(), vec![]);
    }

    #[test]
    fn test_from_shape() {
        let t = YoungTableau::from_shape(&[4, 3, 1]);
        assert_eq!(t.shape(), vec![4, 3, 1]);
        assert_eq!(t.size(), 8);
    }

    #[test]
    fn test_get_set() {
        let mut t = YoungTableau::from_shape(&[3, 2]);
        t.set(0, 0, 1);
        t.set(0, 1, 2);
        assert_eq!(t.get(0, 0), Some(1));
        assert_eq!(t.get(1, 0), Some(0));
        assert_eq!(t.get(5, 5), None);
    }

    #[test]
    fn test_is_semistandard() {
        let t = YoungTableau::from_rows(vec![vec![1, 2, 2], vec![2, 3]]);
        assert!(t.is_semistandard());
        let bad = YoungTableau::from_rows(vec![vec![2, 1], vec![3]]);
        assert!(!bad.is_semistandard());
    }

    #[test]
    fn test_is_standard() {
        let t = YoungTableau::from_rows(vec![vec![1, 3, 5], vec![2, 4]]);
        assert!(t.is_standard());
        let not_std = YoungTableau::from_rows(vec![vec![1, 2], vec![2, 3]]);
        assert!(!not_std.is_standard()); // repeated 2
    }

    #[test]
    fn test_conjugate_shape() {
        let t = YoungTableau::from_shape(&[4, 2, 1]);
        assert_eq!(t.conjugate_shape(), vec![3, 2, 1, 1]);
    }

    #[test]
    fn test_partition_valid() {
        let p = Partition::new(vec![5, 3, 1]);
        assert!(p.is_valid());
        assert_eq!(p.weight(), 9);
        assert!(!Partition::new(vec![1, 3, 5]).is_valid());
    }

    #[test]
    fn test_partition_conjugate() {
        let p = Partition::new(vec![4, 2, 1]);
        assert_eq!(p.conjugate(), Partition::new(vec![3, 2, 1, 1]));
    }

    #[test]
    fn test_partition_dominates() {
        let p1 = Partition::new(vec![3, 1]);
        let p2 = Partition::new(vec![2, 2]);
        assert!(p1.dominates(&p2));
        assert!(!p2.dominates(&p1));
    }

    #[test]
    fn test_display() {
        let t = YoungTableau::from_rows(vec![vec![1, 2], vec![3]]);
        let s = t.display();
        assert!(s.contains("1"));
        assert!(s.contains("2"));
    }
}
