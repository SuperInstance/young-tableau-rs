//! Hook length formula and hook length computations for Young tableaux.

use crate::tableau::Partition;

/// Compute the hook length at cell (row, col) in a tableau of the given shape.
///
/// The hook length is the number of cells to the right in the same row,
/// plus the number of cells below in the same column, plus 1 (the cell itself).
pub fn hook_length(shape: &[usize], row: usize, col: usize) -> u32 {
    if row >= shape.len() || col >= shape[row] {
        return 0;
    }
    let right = (shape[row] - col - 1) as u32;
    let below = shape.iter().skip(row + 1).filter(|&&len| len > col).count() as u32;
    right + below + 1
}

/// Compute all hook lengths for a given shape.
///
/// Returns a 2D vector where `result[row][col]` is the hook length at that cell.
pub fn all_hook_lengths(shape: &[usize]) -> Vec<Vec<u32>> {
    shape
        .iter()
        .enumerate()
        .map(|(row, &len)| (0..len).map(|col| hook_length(shape, row, col)).collect())
        .collect()
}

/// Compute the number of standard Young tableaux of a given shape
/// using the hook length formula:
///
/// `f^λ = n! / ∏ hook(c)`
///
/// where the product is over all cells c of the partition.
pub fn hook_length_count(shape: &[usize]) -> u64 {
    let n: usize = shape.iter().sum();
    if n == 0 {
        return 1;
    }
    let n_fact = factorial(n as u64);
    let hooks_product: u64 = all_hook_lengths(shape)
        .iter()
        .flat_map(|row| row.iter())
        .fold(1u64, |acc, &h| acc * h as u64);

    n_fact / hooks_product
}

/// Compute n!.
fn factorial(n: u64) -> u64 {
    (1..=n).product()
}

/// Compute the hook length formula for a partition.
pub fn syt_count(partition: &Partition) -> u64 {
    hook_length_count(&partition.0)
}

/// Compute the product of all hook lengths.
pub fn hook_product(shape: &[usize]) -> u64 {
    all_hook_lengths(shape)
        .iter()
        .flat_map(|row| row.iter())
        .fold(1u64, |acc, &h| acc * h as u64)
}

/// Compute the hook lengths of the conjugate partition shape.
pub fn conjugate_hook_lengths(shape: &[usize]) -> Vec<Vec<u32>> {
    let conj: Vec<usize> = {
        if shape.is_empty() {
            vec![]
        } else {
            let max_cols = shape[0];
            (0..max_cols).map(|col| shape.iter().filter(|&&len| len > col).count()).collect()
        }
    };
    all_hook_lengths(&conj)
}

/// Compute the dimension of the irreducible representation of S_n
/// corresponding to the given partition. This is the same as syt_count.
pub fn representation_dimension(partition: &Partition) -> u64 {
    syt_count(partition)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hook_length_basic() {
        // Shape [3, 2, 1] (standard hook lengths: [[5,3,1],[3,1],[1]])
        let shape = vec![3, 2, 1];
        assert_eq!(hook_length(&shape, 0, 0), 5);
        assert_eq!(hook_length(&shape, 0, 1), 3);
        assert_eq!(hook_length(&shape, 0, 2), 1);
        assert_eq!(hook_length(&shape, 1, 0), 3);
        assert_eq!(hook_length(&shape, 1, 1), 1);
        assert_eq!(hook_length(&shape, 2, 0), 1);
    }

    #[test]
    fn test_hook_length_out_of_bounds() {
        assert_eq!(hook_length(&[3, 2], 0, 5), 0);
        assert_eq!(hook_length(&[3, 2], 5, 0), 0);
    }

    #[test]
    fn test_all_hook_lengths() {
        let shape = vec![2, 1];
        let hooks = all_hook_lengths(&shape);
        assert_eq!(hooks, vec![vec![3, 1], vec![1]]);
    }

    #[test]
    fn test_hook_length_count_single_row() {
        // Shape [n] has n! / n! = 1 SYT
        assert_eq!(hook_length_count(&[5]), 1);
        assert_eq!(hook_length_count(&[1]), 1);
    }

    #[test]
    fn test_hook_length_count_single_col() {
        // Shape [1,1,1] has 3! / 3! = 1 SYT
        assert_eq!(hook_length_count(&[1, 1, 1]), 1);
    }

    #[test]
    fn test_hook_length_count_2_1() {
        // Shape [2,1]: 3! / (3*1) = 2
        assert_eq!(hook_length_count(&[2, 1]), 2);
    }

    #[test]
    fn test_hook_length_count_2_2() {
        // Shape [2,2]: 4! / (3*2*2*1) = 24/12 = 2
        assert_eq!(hook_length_count(&[2, 2]), 2);
    }

    #[test]
    fn test_hook_length_count_3_2_1() {
        // Shape [3,2,1]: 6! / (5*3*1*3*1*1) = 720/45 = 16
        assert_eq!(hook_length_count(&[3, 2, 1]), 16);
    }

    #[test]
    fn test_hook_product() {
        assert_eq!(hook_product(&[2, 1]), 3);
        assert_eq!(hook_product(&[3, 2, 1]), 45);
    }

    #[test]
    fn test_conjugate_hook_lengths() {
        // Shape [3,1] has conjugate [2,1,1]
        let conj_hooks = conjugate_hook_lengths(&[3, 1]);
        // Conjugate [2,1,1] hook lengths: [[4,1],[2],[1]]
        assert_eq!(conj_hooks, vec![vec![4, 1], vec![2], vec![1]]);
    }

    #[test]
    fn test_representation_dimension() {
        let p = Partition::new(vec![3, 2, 1]);
        assert_eq!(representation_dimension(&p), 16);
    }

    #[test]
    fn test_empty_partition() {
        assert_eq!(hook_length_count(&[]), 1);
        assert_eq!(hook_product(&[]), 1);
    }

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(0), 1);
        assert_eq!(factorial(5), 120);
        assert_eq!(factorial(10), 3628800);
    }
}
