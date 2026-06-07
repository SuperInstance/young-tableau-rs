//! Standard Young tableaux generation and validation.

use crate::hook::hook_length_count;
use crate::tableau::{Partition, YoungTableau};

/// Generate all standard Young tableaux of a given shape using backtracking.
///
/// Returns all SYTs, each containing the numbers 1..=n exactly once.
pub fn generate_standard(shape: &[usize]) -> Vec<YoungTableau> {
    let n: usize = shape.iter().sum();
    if n == 0 {
        return vec![YoungTableau::empty()];
    }
    let mut results = vec![];
    let mut tableau = YoungTableau::from_shape(shape);
    backtrack(&mut tableau, shape, 1, n as u32, &mut results);
    results
}

fn backtrack(
    tableau: &mut YoungTableau,
    shape: &[usize],
    next_val: u32,
    max_val: u32,
    results: &mut Vec<YoungTableau>,
) {
    if next_val > max_val {
        results.push(tableau.clone());
        return;
    }

    for row in 0..shape.len() {
        for col in 0..shape[row] {
            // Cell must be empty (value 0)
            if tableau.get(row, col) != Some(0) {
                continue;
            }

            // Check row constraint: strictly increasing left-to-right
            if col > 0 {
                match tableau.get(row, col - 1) {
                    Some(v) if v == 0 || v >= next_val => continue,
                    _ => {}
                }
            }

            // Check column constraint: strictly increasing top-to-bottom
            if row > 0 {
                match tableau.get(row - 1, col) {
                    Some(v) if v == 0 || v >= next_val => continue,
                    _ => {}
                }
            }

            // For non-first-row cells, the cell above-left must be filled
            // (to ensure canonical filling order)
            if row > 0 && col > 0
                && tableau.get(row - 1, col - 1) == Some(0) {
                    continue;
                }

            tableau.set(row, col, next_val);
            backtrack(tableau, shape, next_val + 1, max_val, results);
            tableau.set(row, col, 0);
        }
    }
}

/// Count the number of standard Young tableaux of the given shape
/// using the hook length formula (verification method).
pub fn count_standard(shape: &[usize]) -> u64 {
    hook_length_count(shape)
}

/// Validate that all generated SYTs are indeed standard.
pub fn validate_all_standard(syts: &[YoungTableau]) -> bool {
    syts.iter().all(|t| t.is_standard())
}

/// Generate all standard Young tableaux of size n across all partitions.
pub fn all_standard_of_size(n: usize) -> Vec<(Partition, Vec<YoungTableau>)> {
    let partitions = generate_partitions(n);
    partitions
        .into_iter()
        .map(|p| {
            let syts = generate_standard(&p.0);
            (p, syts)
        })
        .collect()
}

/// Generate all integer partitions of n.
pub fn generate_partitions(n: usize) -> Vec<Partition> {
    if n == 0 {
        return vec![Partition::new(vec![])];
    }
    let mut result = vec![];
    let mut current = vec![];
    partition_helper(n, n, &mut current, &mut result);
    result
}

fn partition_helper(remaining: usize, max_part: usize, current: &mut Vec<usize>, result: &mut Vec<Partition>) {
    if remaining == 0 {
        result.push(Partition::new(current.clone()));
        return;
    }
    let start = max_part.min(remaining);
    for part in (1..=start).rev() {
        current.push(part);
        partition_helper(remaining - part, part, current, result);
        current.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_standard_shape_1() {
        let syts = generate_standard(&[1]);
        assert_eq!(syts.len(), 1);
        assert!(validate_all_standard(&syts));
    }

    #[test]
    fn test_generate_standard_shape_2_1() {
        let syts = generate_standard(&[2, 1]);
        assert_eq!(syts.len(), 2);
        assert!(validate_all_standard(&syts));
    }

    #[test]
    fn test_generate_standard_shape_2_2() {
        let syts = generate_standard(&[2, 2]);
        assert_eq!(syts.len(), 2);
        assert!(validate_all_standard(&syts));
    }

    #[test]
    fn test_generate_standard_shape_3() {
        let syts = generate_standard(&[3]);
        assert_eq!(syts.len(), 1);
    }

    #[test]
    fn test_generate_standard_shape_3_2_1() {
        let syts = generate_standard(&[3, 2, 1]);
        assert_eq!(syts.len(), 16);
        assert!(validate_all_standard(&syts));
    }

    #[test]
    fn test_count_matches_generation() {
        assert_eq!(generate_standard(&[2, 1]).len() as u64, count_standard(&[2, 1]));
        assert_eq!(generate_standard(&[3, 2, 1]).len() as u64, count_standard(&[3, 2, 1]));
    }

    #[test]
    fn test_generate_partitions() {
        let parts = generate_partitions(4);
        assert_eq!(parts.len(), 5); // [4], [3,1], [2,2], [2,1,1], [1,1,1,1]
    }

    #[test]
    fn test_partitions_sum() {
        let parts = generate_partitions(5);
        for p in &parts {
            assert_eq!(p.weight(), 5);
            assert!(p.is_valid());
        }
    }

    #[test]
    fn test_all_standard_of_size() {
        let all = all_standard_of_size(3);
        let total: usize = all.iter().map(|(_, syts)| syts.len()).sum();
        assert_eq!(total, 4); // n=3 has 4 SYTs total
    }

    #[test]
    fn test_empty_shape() {
        let syts = generate_standard(&[]);
        assert_eq!(syts.len(), 1);
        assert_eq!(syts[0].size(), 0);
    }

    #[test]
    fn test_standard_contain_all_values() {
        let syts = generate_standard(&[3, 2]);
        for t in &syts {
            let mut vals: Vec<u32> = t.rows.iter().flat_map(|r| r.iter().copied()).collect();
            vals.sort();
            assert_eq!(vals, vec![1, 2, 3, 4, 5]);
        }
    }
}
