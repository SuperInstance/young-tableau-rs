//! Schur polynomials.
//!
//! The Schur polynomial `s_λ(x_1, ..., x_k)` is the generating function for
//! semistandard Young tableaux of shape `λ` with entries from `{1,…,k}`:
//!
//! ```text
//! s_λ(x) = Σ_T  x^T
//! ```
//!
//! where `x^T = ∏_{(i,j)∈λ} x_{T(i,j)}`.
//!
//! This module also provides the Jacobi–Trudi identity for verification:
//! `s_λ = det( h_{λ_i - i + j} )` where `h_m` is the complete homogeneous
//! symmetric polynomial.

use crate::tableau::YoungTableau;

/// Evaluate a Schur polynomial `s_λ(x_1,…,x_k)` by enumerating all
/// semistandard Young tableaux of shape `λ` with entries in `1..=k`.
///
/// # Complexity
/// Exponential in the size of `λ`; suitable for small shapes (|λ| ≤ 10).
pub fn schur_polynomial(shape: &[usize], xs: &[f64]) -> f64 {
    if shape.is_empty() {
        return 1.0;
    }
    let k = xs.len();
    let mut total = 0.0;
    let mut tableau = YoungTableau::from_shape(shape);
    backtrack_semistandard(&mut tableau, shape, 0, 0, k, xs, &mut total);
    total
}

/// Backtracking enumeration of semistandard tableaux.
/// Check local constraints after setting a single cell.
/// Zeros represent unfilled cells and are ignored.
fn is_valid_after_set(tableau: &YoungTableau, row: usize, col: usize) -> bool {
    let val = match tableau.get(row, col) {
        Some(v) if v > 0 => v,
        _ => return true,
    };
    // Left neighbor: weakly increasing along rows
    if col > 0 {
        if let Some(left) = tableau.get(row, col - 1) {
            if left > 0 && val < left {
                return false;
            }
        }
    }
    // Above neighbor: strictly increasing down columns
    if row > 0 {
        if let Some(above) = tableau.get(row - 1, col) {
            if above > 0 && val <= above {
                return false;
            }
        }
    }
    // Right neighbor (already filled)
    if let Some(right) = tableau.get(row, col + 1) {
        if right > 0 && val > right {
            return false;
        }
    }
    // Below neighbor (already filled)
    if let Some(below) = tableau.get(row + 1, col) {
        if below > 0 && val >= below {
            return false;
        }
    }
    true
}

fn backtrack_semistandard(
    tableau: &mut YoungTableau,
    shape: &[usize],
    row: usize,
    col: usize,
    max_entry: usize,
    xs: &[f64],
    total: &mut f64,
) {
    if row >= shape.len() {
        *total += monomial(tableau, xs);
        return;
    }

    let next = if col + 1 < shape[row] {
        (row, col + 1)
    } else {
        (row + 1, 0)
    };

    let min_val = if col == 0 {
        1
    } else {
        tableau.get(row, col - 1).unwrap_or(1)
    };

    let min_val = if row == 0 {
        min_val
    } else {
        let above = tableau.get(row - 1, col).unwrap_or(0);
        if above == 0 {
            min_val
        } else {
            min_val.max(above + 1)
        }
    };

    for val in min_val..=max_entry as u32 {
        tableau.set(row, col, val);
        if is_valid_after_set(tableau, row, col) {
            backtrack_semistandard(tableau, shape, next.0, next.1, max_entry, xs, total);
        }
        tableau.set(row, col, 0);
    }
}

/// Compute the monomial weight `x^T` of a filled tableau.
fn monomial(tableau: &YoungTableau, xs: &[f64]) -> f64 {
    let mut prod = 1.0;
    for row in &tableau.rows {
        for &v in row {
            let idx = (v as usize).saturating_sub(1);
            if idx < xs.len() {
                prod *= xs[idx];
            }
        }
    }
    prod
}

/// Evaluate the Schur polynomial at `(1,1,…,1)`.
///
/// This equals the number of semistandard Young tableaux of shape `λ`
/// with entries from `{1,…,k}`.
pub fn schur_at_ones(shape: &[usize], k: usize) -> u64 {
    let xs: Vec<f64> = vec![1.0; k];
    schur_polynomial(shape, &xs) as u64
}

/// The complete homogeneous symmetric polynomial `h_m` evaluated at `xs`.
///
/// `h_m` is the sum of all degree-m monomials in the variables `xs`.
pub fn complete_homogeneous(m: usize, xs: &[f64]) -> f64 {
    if m == 0 {
        return 1.0;
    }
    if xs.is_empty() {
        return 0.0;
    }
    // Dynamic programming: dp[i] = sum of monomials of total degree i
    // Forward iteration because variables can be reused.
    let mut dp = vec![0.0; m + 1];
    dp[0] = 1.0;
    for &x in xs {
        for d in 1..=m {
            dp[d] += x * dp[d - 1];
        }
    }
    dp[m]
}

/// Jacobi–Trudi identity: `s_λ = det( h_{λ_i - i + j} )`.
///
/// Returns the determinant of the `ℓ×ℓ` matrix where `ℓ = len(λ)` and
/// `h_m` is the complete homogeneous symmetric polynomial of degree `m`.
/// For `m < 0` the entry is 0; for `m = 0` it is 1.
pub fn schur_jacobi_trudi(shape: &[usize], xs: &[f64]) -> f64 {
    let n = shape.len();
    if n == 0 {
        return 1.0;
    }
    // Build matrix
    let mut mat = vec![vec![0.0; n]; n];
    for (i, row) in mat.iter_mut().enumerate().take(n) {
        for (j, entry) in row.iter_mut().enumerate().take(n) {
            let m = shape[i] as i64 + j as i64 - i as i64;
            if m >= 0 {
                *entry = complete_homogeneous(m as usize, xs);
            }
        }
    }
    determinant(&mat)
}

/// Determinant of a square matrix via Gaussian elimination.
fn determinant(mat: &[Vec<f64>]) -> f64 {
    let n = mat.len();
    if n == 0 {
        return 1.0;
    }
    let mut a: Vec<Vec<f64>> = mat.iter().map(|r| r.to_vec()).collect();
    let mut det = 1.0;
    let mut sign = 1.0;

    for i in 0..n {
        // Find pivot
        let mut pivot = i;
        for r in i..n {
            if a[r][i].abs() > a[pivot][i].abs() {
                pivot = r;
            }
        }
        if a[pivot][i].abs() < 1e-12 {
            return 0.0;
        }
        if pivot != i {
            a.swap(i, pivot);
            sign = -sign;
        }
        let piv = a[i][i];
        det *= piv;
        for r in i + 1..n {
            let factor = a[r][i] / piv;
            let (above, current) = a.split_at_mut(r);
            let row_i = &above[i];
            for c in i..n {
                current[0][c] -= factor * row_i[c];
            }
        }
    }
    sign * det
}

/// Verify that the combinatorial definition and Jacobi–Trudi agree.
pub fn verify_jacobi_trudi(shape: &[usize], xs: &[f64]) -> bool {
    let s1 = schur_polynomial(shape, xs);
    let s2 = schur_jacobi_trudi(shape, xs);
    (s1 - s2).abs() < 1e-6
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_schur_empty_shape() {
        assert_eq!(schur_polynomial(&[], &[1.0, 2.0]), 1.0);
    }

    #[test]
    fn test_schur_single_row() {
        // s_{(m)} = h_m
        let xs = vec![1.0, 2.0, 3.0];
        assert!((schur_polynomial(&[2], &xs) - complete_homogeneous(2, &xs)).abs() < 1e-10);
    }

    #[test]
    fn test_schur_single_column() {
        // s_{(1^m)} = e_m (elementary symmetric), but for m=2:
        // e_2(x1,x2,x3) = x1x2 + x1x3 + x2x3
        let xs = vec![1.0, 2.0, 3.0];
        let s = schur_polynomial(&[1, 1], &xs);
        let expected = 1.0 * 2.0 + 1.0 * 3.0 + 2.0 * 3.0;
        assert!((s - expected).abs() < 1e-10);
    }

    #[test]
    fn test_schur_at_ones_2_1() {
        // Number of SSYT of shape (2,1) with entries from {1,2,3}
        // = 8
        assert_eq!(schur_at_ones(&[2, 1], 3), 8);
    }

    #[test]
    fn test_jacobi_trudi_2_1() {
        let shape = vec![2, 1];
        let xs = vec![1.0, 2.0, 3.0];
        assert!(verify_jacobi_trudi(&shape, &xs));
    }

    #[test]
    fn test_jacobi_trudi_3_2() {
        let shape = vec![3, 2];
        let xs = vec![0.5, 1.5, 2.0, 0.3];
        assert!(verify_jacobi_trudi(&shape, &xs));
    }

    #[test]
    fn test_complete_homogeneous() {
        let xs = vec![1.0, 2.0];
        assert!((complete_homogeneous(0, &xs) - 1.0).abs() < 1e-10);
        assert!((complete_homogeneous(1, &xs) - 3.0).abs() < 1e-10);
        assert!((complete_homogeneous(2, &xs) - (1.0 + 2.0 + 4.0)).abs() < 1e-10);
    }

    #[test]
    fn test_determinant_identity() {
        let mat = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        assert!((determinant(&mat) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_determinant_2x2() {
        let mat = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        assert!((determinant(&mat) - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_schur_at_ones_single_row() {
        // s_{(3)}(1,1,1) = h_3(1,1,1) = number of ways to write 3 as ordered sum of 3 vars
        // = 10 (stars and bars)
        assert_eq!(schur_at_ones(&[3], 3), 10);
    }

    #[test]
    fn test_schur_homogeneous_consistency() {
        // For a single-row shape (m), s_{(m)} = h_m
        for m in 1..=4 {
            let xs = vec![1.0, 2.0];
            let s = schur_polynomial(&[m], &xs);
            let h = complete_homogeneous(m, &xs);
            assert!((s - h).abs() < 1e-8, "Mismatch for m={}: {} vs {}", m, s, h);
        }
    }
}
