//! # young-tableau-rs
//!
//! A library for working with Young tableaux: Robinson-Schensted insertion,
//! bumping, hook length formula, and Schensted correspondence.
//!
//! ## Modules
//!
//! - [`tableau`] — Core Young tableau data structure
//! - [`insertion`] — Robinson-Schensted row insertion and bumping
//! - [`hook`] — Hook length formula and hook lengths
//! - [`standard`] — Standard Young tableaux generation and validation
//! - [`correspondence`] — Robinson-Schensted-Knuth correspondence

pub mod correspondence;
pub mod hook;
pub mod insertion;
pub mod schur;
pub mod standard;
pub mod tableau;

pub use tableau::YoungTableau;
