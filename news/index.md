# Changelog

## survdistr 0.0.2

- Added
  [`convert_to_surv()`](https://survdistr.mlr-org.com/reference/convert_to_surv.md)
  for efficient transformation of hazard/density functions to survival
  at anchor time points (uses C++ internally).
- Introduced
  [`trim_duplicates()`](https://survdistr.mlr-org.com/reference/trim_duplicates.md)
  to remove repeated S(t) values across the time axis, with tolerance.
- Added
  [`interp_cif()`](https://survdistr.mlr-org.com/reference/interp_cif.md)
  for constant interpolation of CIF (uses C++ internally).
- Added
  [`extract_times()`](https://survdistr.mlr-org.com/reference/extract_times.md)
  to consistently obtain and validate time points from vectors or
  matrices.
- Unified assertion function to one:
  [`assert_prob()`](https://survdistr.mlr-org.com/reference/assert_prob.md).
- Unified `mat_interp()` and `vec_interp()` into a single
  [`interp()`](https://survdistr.mlr-org.com/reference/interp.md)
  function.
- Enhanced
  [`interp()`](https://survdistr.mlr-org.com/reference/interp.md):
  - Added aliases for the three interpolation options (e.g., `const_haz`
    is equivalent to `exp_surv`).
  - Implemented S(t) interpolation in C++ for all 3 options
    (`output = "surv"`); H(t) and F(t) follow similarly.
  - Implemented f(t) interpolation in C++ for all 3 options
    (`output = "density"`)
  - Implemented h(t) interpolation in C++ for all 3 options
    (`output = "hazard"`)
  - Requires unique and ordered `eval_times`, resulting in significant
    C++ speed-up.
  - Allows passing `times` separately to preserve precision.
- Improved `survDistr`:
  - Removed `data_type` argument; input now only S(t).
  - Added aliases for the three interpolation options.
  - Allows passing `times` during construction separately for precision.
  - Introduced [`filter()`](https://rdrr.io/r/stats/filter.html) method
    for in-place object filtering.
  - Added `rows` argument to filter input matrix before interpolation
    without altering the object.
  - `method` and `times` are now active fields, not public members.
  - Added S3 conversion via `as_survDistr(x)`.
- Added unit tests for all features.

## survdistr 0.0.1

- Base `survDistr` class.
- `assert_prob_matrix()` for asserting survival, CDF and CIF matrices.
- `mat_interp()` and `vec_interp()` for constant and linear
  interpolation of survival, CDF and CIF matrices/vectors (using
  internal `Rcpp` methods).
