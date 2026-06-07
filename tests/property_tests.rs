use proptest::prelude::*;
use young_tableau_rs::correspondence::{rsk, rsk_inverse};
use young_tableau_rs::hook::{hook_length_count, syt_count};
use young_tableau_rs::insertion::{insertion_tableau, recording_tableau};
use young_tableau_rs::schur::{schur_jacobi_trudi, schur_polynomial, verify_jacobi_trudi};
use young_tableau_rs::tableau::Partition;

prop_compose! {
    fn vec_permutation()(n in 1usize..=8) -> Vec<u32> {
        let mut v: Vec<u32> = (1..=n as u32).collect();
        let mut seed = n as u64 * 123456789;
        // simple xorshift shuffle
        for i in (1..v.len()).rev() {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            let j = (seed % (i as u64 + 1)) as usize;
            v.swap(i, j);
        }
        v
    }
}

prop_compose! {
    fn small_partition()(mut parts in prop::collection::vec(1usize..=5, 0..=5)) -> Vec<usize> {
        parts.sort_by(|a, b| b.cmp(a));
        parts
    }
}

prop_compose! {
    fn small_vars()(n in 1usize..=4) -> Vec<f64> {
        (0..n).map(|i| 0.1 + 0.5 * (i as f64 + 1.0)).collect()
    }
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 128,
        ..ProptestConfig::default()
    })]

    #[test]
    fn prop_rsk_roundtrip(perm in vec_permutation()) {
        let pair = rsk(&perm);
        let recovered = rsk_inverse(&pair);
        prop_assert_eq!(perm, recovered);
    }

    #[test]
    fn prop_insertion_recording_same_shape(perm in vec_permutation()) {
        let p = insertion_tableau(&perm);
        let q = recording_tableau(&perm);
        prop_assert_eq!(p.shape(), q.shape());
    }

    #[test]
    fn prop_hook_length_count_matches_syt(shape in small_partition()) {
        let n: usize = shape.iter().sum();
        if n == 0 {
            prop_assert_eq!(hook_length_count(&shape), 1);
            return Ok(());
        }
        let count = hook_length_count(&shape);
        let partition = Partition::new(shape.clone());
        prop_assert_eq!(count, syt_count(&partition));
    }

    #[test]
    fn prop_schur_jacobi_trudi_agreement(shape in small_partition(), xs in small_vars()) {
        if shape.is_empty() {
            prop_assert!(verify_jacobi_trudi(&shape, &xs));
            return Ok(());
        }
        let n: usize = shape.iter().sum();
        // skip shapes too large for enumeration (limit to |λ| ≤ 6)
        if n > 6 {
            return Ok(());
        }
        prop_assert!(verify_jacobi_trudi(&shape, &xs));
        let s1 = schur_polynomial(&shape, &xs);
        let s2 = schur_jacobi_trudi(&shape, &xs);
        prop_assert!((s1 - s2).abs() < 1e-6, "Jacobi–Trudi mismatch: {} vs {}", s1, s2);
    }

    #[test]
    fn prop_schur_homogeneity(shape in small_partition(), xs in small_vars()) {
        let n: usize = shape.iter().sum();
        if n == 0 || n > 6 {
            return Ok(());
        }
        let s = schur_polynomial(&shape, &xs);
        let t: f64 = 2.0;
        let scaled_xs: Vec<f64> = xs.iter().map(|&x| t * x).collect();
        let scaled_s = schur_polynomial(&shape, &scaled_xs);
        let expected = t.powi(n as i32) * s;
        prop_assert!((scaled_s - expected).abs() < 1e-6 * expected.abs().max(1.0),
            "Homogeneity failed: {} vs {}", scaled_s, expected);
    }
}
