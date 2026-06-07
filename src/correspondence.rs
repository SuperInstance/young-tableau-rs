//! Robinson-Schensted-Knuth correspondence: bijection between permutations
//! and pairs of standard Young tableaux (insertion and recording).

use crate::insertion::{insertion_tableau, recording_tableau};
use crate::tableau::YoungTableau;

/// A pair of tableaux from the RSK correspondence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RSKPair {
    /// The insertion (P) tableau.
    pub p: YoungTableau,
    /// The recording (Q) tableau.
    pub q: YoungTableau,
}

/// Compute the RSK correspondence for a permutation.
///
/// Given a permutation π of [1..n], returns the pair (P, Q) of standard Young tableaux
/// where P is the insertion tableau and Q is the recording tableau.
pub fn rsk(perm: &[u32]) -> RSKPair {
    assert!(!perm.is_empty(), "permutation must be non-empty");
    let p = insertion_tableau(perm);
    let q = recording_tableau(perm);
    RSKPair { p, q }
}

/// Recover the permutation from an RSK pair by reverse-insertion.
///
/// This is the inverse of the RSK correspondence.
pub fn rsk_inverse(pair: &RSKPair) -> Vec<u32> {
    let n = pair.p.size();
    if n == 0 {
        return vec![];
    }

    let mut p = pair.p.clone();
    let q = pair.q.clone();
    let mut result = vec![0u32; n];

    // Read values from Q to determine reverse order
    let mut q_values: Vec<(u32, usize, usize)> = vec![];
    for (r, row) in q.rows.iter().enumerate() {
        for (c, &v) in row.iter().enumerate() {
            q_values.push((v, r, c));
        }
    }
    q_values.sort_by(|a, b| b.0.cmp(&a.0)); // Descending order

    for &(val, _, _) in &q_values {
        let idx = (val - 1) as usize;
        // Find a corner cell of P
        let (cr, cc) = find_corner(&p);
        let (new_p, bumped) = crate::insertion::reverse_insert(&p, cr, cc);
        result[idx] = bumped;
        p = new_p;
    }

    result
}

/// Find a corner cell (maximal in both coordinates) of the tableau.
fn find_corner(tableau: &YoungTableau) -> (usize, usize) {
    let last_row = tableau.num_rows() - 1;
    let last_col = tableau.rows[last_row].len() - 1;
    (last_row, last_col)
}

/// Compute the length of the longest increasing subsequence of a permutation.
///
/// This equals the number of columns in the insertion tableau.
pub fn longest_increasing_subsequence_length(perm: &[u32]) -> usize {
    if perm.is_empty() {
        return 0;
    }
    let p = insertion_tableau(perm);
    p.rows[0].len()
}

/// Compute the length of the longest decreasing subsequence of a permutation.
///
/// This equals the number of rows in the insertion tableau.
pub fn longest_decreasing_subsequence_length(perm: &[u32]) -> usize {
    if perm.is_empty() {
        return 0;
    }
    let p = insertion_tableau(perm);
    p.num_rows()
}

/// Compute the insertion tableau shape for a permutation.
pub fn insertion_shape(perm: &[u32]) -> Vec<usize> {
    insertion_tableau(perm).shape()
}

/// Check if two RSK pairs are consistent (P and Q have the same shape).
pub fn is_consistent(pair: &RSKPair) -> bool {
    pair.p.shape() == pair.q.shape()
}

/// Compute the RSK correspondence for a general sequence (not necessarily a permutation).
pub fn rsk_general(sequence: &[u32]) -> RSKPair {
    let p = insertion_tableau(sequence);
    let q = recording_tableau(sequence);
    RSKPair { p, q }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rsk_identity() {
        let pair = rsk(&[1, 2, 3, 4]);
        assert!(pair.p.is_standard());
        assert!(pair.q.is_standard());
        assert!(is_consistent(&pair));
    }

    #[test]
    fn test_rsk_reverse() {
        let pair = rsk(&[4, 3, 2, 1]);
        assert!(pair.p.is_standard());
        assert!(pair.q.is_standard());
        // Reverse permutation should give a single-column tableau
        assert_eq!(pair.p.num_rows(), 4);
    }

    #[test]
    fn test_rsk_consistency() {
        let pair = rsk(&[3, 1, 4, 2]);
        assert!(is_consistent(&pair));
    }

    #[test]
    fn test_lis_length() {
        assert_eq!(longest_increasing_subsequence_length(&[1, 2, 3]), 3);
        assert_eq!(longest_increasing_subsequence_length(&[3, 2, 1]), 1);
        assert_eq!(longest_increasing_subsequence_length(&[2, 1, 3]), 2);
    }

    #[test]
    fn test_lds_length() {
        assert_eq!(longest_decreasing_subsequence_length(&[3, 2, 1]), 3);
        assert_eq!(longest_decreasing_subsequence_length(&[1, 2, 3]), 1);
    }

    #[test]
    fn test_lis_lds_duality() {
        let perm = vec![3, 1, 4, 1, 5, 9, 2, 6];
        let lis = longest_increasing_subsequence_length(&perm);
        let lds = longest_decreasing_subsequence_length(&perm);
        assert!(lis >= 1);
        assert!(lds >= 1);
    }

    #[test]
    fn test_insertion_shape() {
        let shape = insertion_shape(&[3, 1, 2]);
        assert!(!shape.is_empty());
        let total: usize = shape.iter().sum();
        assert_eq!(total, 3);
    }

    #[test]
    fn test_rsk_general() {
        let pair = rsk_general(&[1, 2, 1, 3, 2]);
        assert!(is_consistent(&pair));
        assert_eq!(pair.p.size(), 5);
    }

    #[test]
    fn test_empty_permutation() {
        assert_eq!(longest_increasing_subsequence_length(&[]), 0);
        assert_eq!(longest_decreasing_subsequence_length(&[]), 0);
    }

    #[test]
    fn test_rsk_pair_shapes() {
        for perm in [[1u32, 2, 3], [3, 2, 1], [2, 1, 3], [1, 3, 2]] {
            let pair = rsk(&perm);
            assert_eq!(pair.p.shape(), pair.q.shape(), "Shape mismatch for {:?}", perm);
        }
    }

    #[test]
    fn test_corner_finding() {
        let t = YoungTableau::from_rows(vec![vec![1, 2, 3], vec![4, 5]]);
        let (r, c) = find_corner(&t);
        assert_eq!(r, 1);
        assert_eq!(c, 1);
    }
}
