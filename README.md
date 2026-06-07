# young-tableau-rs

Young tableaux operations: Robinson-Schensted insertion, bumping, hook length formula, and Schensted correspondence.

## Features

- **Tableau**: Core Young tableau data structure with semistandard/standard validation
- **Insertion**: Robinson-Schensted row insertion and bumping paths
- **Hook Length**: Hook length formula for counting standard Young tableaux
- **Standard**: Generation of all SYTs of a given shape
- **Correspondence**: RSK correspondence, LIS/LDS computation

Pure Rust, no external dependencies.

## Usage

```rust
use young_tableau_rs::insertion::insertion_tableau;

let t = insertion_tableau(&[3, 1, 4, 1, 5]);
assert_eq!(t.size(), 5);
assert!(t.is_semistandard());
```

License: MIT OR Apache-2.0
