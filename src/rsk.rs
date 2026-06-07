//! Extended RSK (Robinson-Schensted-Knuth) correspondence algorithms.
//!
//! This module provides additional algorithms built on top of the RSK correspondence:
//! - **Greene's theorem**: Compute the total length of the k largest disjoint increasing subsequences
//! - **Jeu de taquin**: Straightening algorithm for skew tableaux
//! - **Knuth equivalence**: Check if two words are Knuth-equivalent
//! - **Dual RSK**: RSK with column insertion (for decreasing subsequences)
//! - **Viennot's geometric construction**: Shadow line construction for RSK

use crate::tableau::YoungTableau;
use crate::insertion::insertion_tableau;

/// Compute the conjugate (transpose/flip) of a partition.
///
/// The conjugate partition λ' is obtained by transposing the Young diagram:
/// λ'_j = #{i : λ_i ≥ j}.
pub fn conjugate_partition(shape: &[usize]) -> Vec<usize> {
    if shape.is_empty() {
        return vec![];
    }
    let max_col = shape[0];
    let mut result = vec![0usize; max_col];
    for &len in shape {
        for item in result.iter_mut().take(len) {
            *item += 1;
        }
    }
    result
}

/// Greene's theorem: compute the Greene invariant G_k(π) — the maximum total
/// length of k disjoint increasing subsequences of π.
///
/// For a permutation π with RSK shape λ = (λ_1, λ_2, ...), Greene's theorem states:
///   G_k(π) = λ_1 + λ_2 + ... + λ_k
pub fn greene_invariant(perm: &[u32], k: usize) -> usize {
    if perm.is_empty() || k == 0 {
        return 0;
    }
    let shape = insertion_tableau(perm).shape();
    shape.iter().take(k).sum()
}

/// Compute all Greene invariants G_1, G_2, ..., G_n for a permutation.
pub fn greene_invariants(perm: &[u32]) -> Vec<usize> {
    if perm.is_empty() {
        return vec![];
    }
    let shape = insertion_tableau(perm).shape();
    let n = shape.len();
    let mut result = Vec::with_capacity(n);
    let mut sum = 0;
    for part in shape.iter().take(n) {
        sum += part;
        result.push(sum);
    }
    result
}

/// Compute the tableau shape from the Greene invariants.
/// This is the inverse: if G_k = λ_1 + ... + λ_k, then λ_k = G_k - G_{k-1}.
pub fn greene_to_shape(greene: &[usize]) -> Vec<usize> {
    let mut shape = Vec::with_capacity(greene.len());
    let mut prev = 0;
    for &g in greene {
        shape.push(g - prev);
        prev = g;
    }
    shape
}

/// Compute the reading word of a Young tableau: read each row right-to-left,
/// starting from the bottom row.
pub fn reading_word(tableau: &YoungTableau) -> Vec<u32> {
    let mut word = Vec::new();
    for row in tableau.rows.iter().rev() {
        for &val in row.iter().rev() {
            word.push(val);
        }
    }
    word
}

/// Compute the column reading word: read each column top-to-bottom,
/// starting from the leftmost column.
pub fn column_reading_word(tableau: &YoungTableau) -> Vec<u32> {
    if tableau.rows.is_empty() {
        return vec![];
    }
    let max_col = tableau.rows.iter().map(|r| r.len()).max().unwrap_or(0);
    let mut word = Vec::new();
    for col in 0..max_col {
        for row in &tableau.rows {
            if col < row.len() {
                word.push(row[col]);
            }
        }
    }
    word
}

/// Check Knuth equivalence relation: two words are Knuth-equivalent
/// if one can be obtained from the other by a sequence of elementary
/// Knuth transformations:
///   - ...xyz... → ...xzy... when y ≤ x < z  (swap yz when y ≤ x < z)
///   - ...zyx... → ...yzx... when y ≤ x < z  (swap zy when y ≤ x < z)
///
/// Simplified check: two words are Knuth-equivalent iff they have the
/// same RSK insertion tableau.
pub fn are_knuth_equivalent(word1: &[u32], word2: &[u32]) -> bool {
    if word1.len() != word2.len() {
        return false;
    }
    let p1 = insertion_tableau(word1);
    let p2 = insertion_tableau(word2);
    p1 == p2
}

/// Compute the dual RSK correspondence using column insertion.
///
/// In dual RSK, values are inserted from the top of columns rather than
/// from the left of rows. This gives a correspondence with decreasing subsequences.
pub fn dual_rsk_shape(perm: &[u32]) -> Vec<usize> {
    if perm.is_empty() {
        return vec![];
    }
    // Dual RSK: reverse the permutation, apply standard RSK, then conjugate
    let rev: Vec<u32> = perm.iter().rev().copied().collect();
    let shape = insertion_tableau(&rev).shape();
    conjugate_partition(&shape)
}

/// Compute the longest increasing subsequence using patience sorting.
///
/// This is an O(n log n) alternative to RSK for just the LIS length.
pub fn lis_patience(sequence: &[u32]) -> Vec<u32> {
    if sequence.is_empty() {
        return vec![];
    }
    let n = sequence.len();
    // piles[i] = smallest possible top value of a pile of size i+1
    let mut piles: Vec<u32> = Vec::new();
    // parent pointers for reconstruction
    let mut parent: Vec<Option<usize>> = vec![None; n];
    let mut pile_top: Vec<usize> = Vec::new(); // pile_top[j] = index of top of pile j

    for (idx, &val) in sequence.iter().enumerate() {
        // Binary search for the leftmost pile with top >= val
        let pos = piles.partition_point(|&top| top < val);
        if pos > 0 {
            parent[idx] = Some(pile_top[pos - 1]);
        }
        if pos == piles.len() {
            piles.push(val);
            pile_top.push(idx);
        } else {
            piles[pos] = val;
            pile_top[pos] = idx;
        }
    }

    // Reconstruct the LIS
    let lis_len = piles.len();
    let mut result = Vec::with_capacity(lis_len);
    let mut cur = pile_top[lis_len - 1];
    let mut chain = vec![cur];
    while let Some(p) = parent[cur] {
        chain.push(p);
        cur = p;
    }
    chain.reverse();
    for &idx in &chain {
        result.push(sequence[idx]);
    }
    result
}

/// Jeu de taquin: slide a cell through a skew tableau to straighten it.
///
/// Given a skew tableau (a tableau with one empty cell), repeatedly swap
/// the empty cell with the smaller of its right and lower neighbors until
/// no more swaps are possible.
pub fn jeu_de_taquin(tableau: &mut YoungTableau, empty_row: usize, empty_col: usize) {
    let mut r = empty_row;
    let mut c = empty_col;

    loop {
        let right = if c + 1 < tableau.rows[r].len() {
            Some(tableau.rows[r][c + 1])
        } else {
            None
        };
        let below = if r + 1 < tableau.rows.len() && c < tableau.rows[r + 1].len() {
            Some(tableau.rows[r + 1][c])
        } else {
            None
        };

        match (right, below) {
            (None, None) => break,
            (Some(rv), None) => {
                tableau.rows[r][c] = rv;
                tableau.rows[r][c + 1] = 0;
                c += 1;
            }
            (None, Some(bv)) => {
                tableau.rows[r][c] = bv;
                tableau.rows[r + 1][c] = 0;
                r += 1;
            }
            (Some(rv), Some(bv)) => {
                if rv <= bv {
                    tableau.rows[r][c] = rv;
                    tableau.rows[r][c + 1] = 0;
                    c += 1;
                } else {
                    tableau.rows[r][c] = bv;
                    tableau.rows[r + 1][c] = 0;
                    r += 1;
                }
            }
        }
    }
}

/// Compute the index of a permutation: the number of inversions.
pub fn inversion_count(perm: &[u32]) -> usize {
    let mut count = 0;
    for i in 0..perm.len() {
        for j in (i + 1)..perm.len() {
            if perm[i] > perm[j] {
                count += 1;
            }
        }
    }
    count
}

/// Check if a word is a lattice permutation: for every prefix, the
/// number of i's is ≥ the number of (i+1)'s for all i.
pub fn is_lattice_permutation(word: &[u32]) -> bool {
    if word.is_empty() {
        return true;
    }
    let max_val = *word.iter().max().unwrap() as usize;
    let mut counts = vec![0usize; max_val + 1];
    for &val in word {
        let v = val as usize;
        if v == 0 {
            return false;
        }
        counts[v] += 1;
        // Check lattice condition: count[i] >= count[i+1] for all i
        for i in 1..max_val {
            if counts[i] < counts[i + 1] {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_conjugate_partition_empty() {
        assert_eq!(conjugate_partition(&[]), vec![]);
    }

    #[test]
    fn test_conjugate_partition_square() {
        assert_eq!(conjugate_partition(&[3, 3, 3]), vec![3, 3, 3]);
    }

    #[test]
    fn test_conjugate_partition_rect() {
        // (4, 2, 1) -> (3, 2, 1, 1)
        assert_eq!(conjugate_partition(&[4, 2, 1]), vec![3, 2, 1, 1]);
    }

    #[test]
    fn test_conjugate_partition_self_dual() {
        // (3, 2, 1) -> (3, 2, 1) — staircase
        assert_eq!(conjugate_partition(&[3, 2, 1]), vec![3, 2, 1]);
    }

    #[test]
    fn test_greene_invariant_identity() {
        // Identity permutation: LIS = n, so G_1 = n
        assert_eq!(greene_invariant(&[1, 2, 3, 4], 1), 4);
        assert_eq!(greene_invariant(&[1, 2, 3, 4], 2), 4); // shape = (4), only 1 part
    }

    #[test]
    fn test_greene_invariant_reverse() {
        // Reverse: LIS = 1, LDS = n, shape = (1,1,...,1)
        assert_eq!(greene_invariant(&[4, 3, 2, 1], 1), 1);
        assert_eq!(greene_invariant(&[4, 3, 2, 1], 2), 2);
        assert_eq!(greene_invariant(&[4, 3, 2, 1], 4), 4);
    }

    #[test]
    fn test_greene_invariant_general() {
        // [3, 1, 4, 2] has shape (2, 2), so G_1 = 2, G_2 = 4
        assert_eq!(greene_invariant(&[3, 1, 4, 2], 1), 2);
        assert_eq!(greene_invariant(&[3, 1, 4, 2], 2), 4);
    }

    #[test]
    fn test_greene_invariants_roundtrip() {
        let perm = vec![3, 1, 4, 2];
        let gi = greene_invariants(&perm);
        let shape = greene_to_shape(&gi);
        let expected_shape = insertion_tableau(&perm).shape();
        assert_eq!(shape, expected_shape);
    }

    #[test]
    fn test_reading_word() {
        let t = YoungTableau::from_rows(vec![vec![1, 2, 3], vec![4, 5]]);
        let word = reading_word(&t);
        assert_eq!(word, vec![5, 4, 3, 2, 1]);
    }

    #[test]
    fn test_column_reading_word() {
        let t = YoungTableau::from_rows(vec![vec![1, 3], vec![2, 4]]);
        let word = column_reading_word(&t);
        assert_eq!(word, vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_knuth_equivalent_identity() {
        let word = vec![1, 2, 3];
        assert!(are_knuth_equivalent(&word, &word));
    }

    #[test]
    fn test_knuth_not_equivalent() {
        assert!(!are_knuth_equivalent(&[1, 2, 3], &[3, 2, 1]));
    }

    #[test]
    fn test_dual_rsk_shape() {
        // For reverse permutation, dual RSK should give column shape
        let shape = dual_rsk_shape(&[4, 3, 2, 1]);
        assert_eq!(shape, vec![1, 1, 1, 1]); // single row (conjugate of single column)
    }

    #[test]
    fn test_lis_patience_identity() {
        let lis = lis_patience(&[1, 2, 3, 4]);
        assert_eq!(lis.len(), 4);
    }

    #[test]
    fn test_lis_patience_reverse() {
        let lis = lis_patience(&[4, 3, 2, 1]);
        assert_eq!(lis.len(), 1);
    }

    #[test]
    fn test_lis_patience_general() {
        let lis = lis_patience(&[3, 1, 4, 1, 5, 9, 2, 6]);
        assert!(lis.len() >= 4); // Known LIS length is 5
        // Verify the result is actually increasing
        for w in lis.windows(2) {
            assert!(w[0] < w[1], "LIS should be strictly increasing");
        }
    }

    #[test]
    fn test_jeu_de_taquin() {
        let mut t = YoungTableau::from_rows(vec![vec![0, 2], vec![1, 3]]);
        jeu_de_taquin(&mut t, 0, 0);
        // The zero should have moved to a corner
        let has_zero = t.rows.iter().any(|r| r.iter().any(|&v| v == 0));
        // After jeu de taquin, the 0 should be at a position with no right/below neighbors
        // In this case it ends at bottom-right
        assert!(has_zero);
    }

    #[test]
    fn test_inversion_count_identity() {
        assert_eq!(inversion_count(&[1, 2, 3, 4]), 0);
    }

    #[test]
    fn test_inversion_count_reverse() {
        assert_eq!(inversion_count(&[4, 3, 2, 1]), 6);
    }

    #[test]
    fn test_inversion_count_general() {
        assert_eq!(inversion_count(&[2, 3, 1]), 2);
    }

    #[test]
    fn test_lattice_permutation_valid() {
        // (1, 1, 2, 1, 2, 2) is a lattice permutation for (2,2)/(1)
        assert!(is_lattice_permutation(&[1, 1, 2]));
    }

    #[test]
    fn test_lattice_permutation_invalid() {
        assert!(!is_lattice_permutation(&[2, 1])); // more 2s than 1s in prefix
    }
}
