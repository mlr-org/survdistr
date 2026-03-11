# survdistr 0.0.2

- Added `convert_to_surv()` for efficient transformation of hazard/density functions to survival at anchor time points (uses C++ internally).
- Introduced `trim_duplicates()` to remove repeated S(t) values across the time axis, with tolerance.
- Added `interp_cif()` for constant interpolation of CIF (uses C++ internally).
- Unified assertion function to one: `assert_prob()`.
- Unified `mat_interp()` and `vec_interp()` into a single `interp()` function.
- Enhanced `interp()`:
  - Added aliases for the three interpolation options (e.g., `const_haz` is equivalent to `exp_surv`).
  - Implemented S(t) interpolation in C++ for all options; H(t) and F(t) follow similarly.
  - Requires unique and ordered `eval_times`, resulting in significant C++ speed-up.
  - Allows passing `times` separately to preserve precision.
- Improved `survDistr`:
  - Removed `data_type` argument; input now only S(t).
  - Added aliases for the three interpolation options.
  - Allows passing `times` during construction separately for precision.
  - Introduced `filter()` method for in-place object filtering.
  - Added `rows` argument to filter input matrix before interpolation without altering the object.
  - `method` and `times` are now active fields, not public members.
  - Added S3 conversion via `as_survDistr(x)`.
- Added unit tests for all features.

# survdistr 0.0.1

- Base `survDistr` class.
- `assert_prob_matrix()` for asserting survival, CDF and CIF matrices.
- `mat_interp()` and `vec_interp()` for constant and linear interpolation of survival, CDF and CIF matrices/vectors (using internal `Rcpp` methods).
