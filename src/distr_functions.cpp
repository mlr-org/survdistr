#include <Rcpp.h>
using namespace Rcpp;

// Row value interpolation of a survival/cdf matrix
// * `x` is a survival (or cdf) matrix
// * `times` must be non-negative, unique, increasing numbers
// * `new_times` can be unordered and duplicated numeric values
// TODO: new_times ordered will speed up things, can do column ordering in R
// TODO: handle here if new_times = NULL (simply return x)
// * `surv` denotes the type of `x` (survival matrix by default)
// * `inter_type` denotes the interpolation type (default: constant, otherwise linear)
// [[Rcpp::export]]
NumericMatrix rcpp_mat_interp(const NumericMatrix& x,
                              const NumericVector& times,
                              const NumericVector& new_times,
                              bool surv = true,
                              bool constant = true) {
  // Constants for default values
  constexpr double SURV_DEFAULT = 1.0;
  constexpr double CDF_DEFAULT = 0.0;

  // x => [obs x times], S(t) or CDF(t)
  int n_rows = x.nrow(); // Number of observations
  int n_times = x.ncol(); // Number of original time points
  int n_times_new = new_times.length(); // Number of new time points

  // NOTE: `times` is the colnames of `x` matrix
  if (n_times != times.length()) {
    stop("Number of columns in the input matrix must match the length of 'times'.");
  }

  NumericMatrix mat(n_rows, n_times_new); // Resulting matrix: [n_rows x n_times_new]

  for (int i = 0; i < n_rows; i++) { // Iterate over observations
    for (int k = 0; k < n_times_new; k++) { // Iterate over new time points
      double t_new = new_times[k];

      // new time is before the first time point (extrapolation)
      if (t_new < times[0]) {
        if (constant) {
          // Constant interpolation: use default value
          mat(i, k) = surv ? SURV_DEFAULT : CDF_DEFAULT;
        } else {
          // Linear extrapolation considering that S(t = 0) = 1
          double x1 = surv ? SURV_DEFAULT : CDF_DEFAULT;
          double x2 = x(i, 0);

          mat(i, k) = x1 + t_new * (x2 - x1) / times[0];
        }
        continue;
      }

      // new time is after the last time point (extrapolation)
      if (t_new >= times[n_times - 1]) {
        if (constant || n_times <= 1) {
          mat(i, k) = x(i, n_times - 1); // Use last value for constant extrapolation
        } else {
          // Linear extrapolation using the last time interval
          double t1 = times[n_times - 2], t2 = times[n_times - 1];
          double x1 = x(i, n_times - 2), x2 = x(i, n_times - 1);
          double extrapolated_value = x2 + (t_new - t2) * (x2 - x1) / (t2 - t1);

          if (surv) { // Adjust bounds for survival function
            mat(i, k) = std::max(0.0, extrapolated_value);
          } else { // Adjust bounds for CDF
            mat(i, k) = std::min(1.0, extrapolated_value);
          }
        }
        continue;
      }

      // Find the correct interval for interpolation
      for (int j = 0; j < n_times - 1; j++) {
        if (t_new >= times[j] && t_new < times[j + 1]) {
          if (constant) {
            // Constant interpolation (use value at the lower bound)
            mat(i, k) = x(i, j);
          } else {
            // Linear interpolation
            double t1 = times[j], t2 = times[j + 1];
            double x1 = x(i, j), x2 = x(i, j + 1);
            mat(i, k) = x1 + (t_new - t1) * (x2 - x1) / (t2 - t1);
          }

          break;
        }
      }
    }
  }

  return mat;
}

// This function does row-wise (time-)weighted cumulative sum
// NOTE: `times` must be non-negative, unique, increasing numbers
// [[Rcpp::export]]
NumericMatrix rcpp_wcumsum_mat(const NumericMatrix& x,
                               const NumericVector& times) {
  int n_rows = x.nrow();
  int n_times = x.ncol();
  NumericMatrix mat(n_rows, n_times);

  // NOTE: `times` is the colnames of `x` matrix
  if (n_times != times.length()) {
    stop("Number of columns in the input matrix must match the length of 'times'.");
  }

  // Compute time differences (weights)
  std::vector<double> delta_t(n_times);
  delta_t[0] = times[0]; // First interval from time 0
  for (int j = 1; j < n_times; j++) {
    delta_t[j] = times[j] - times[j - 1];
  }

  // Compute cumulative weighted sums row-wise
  for (int i = 0; i < n_rows; i++) {
    mat(i, 0) = x(i, 0) * delta_t[0]; // First column
    for (int j = 1; j < n_times; j++) {
      mat(i, j) = mat(i, j - 1) + x(i, j) * delta_t[j];
    }
  }

  return mat;
}

// This function does row-wise (time-)weighted cumulative sum for in-between time points
// NOTE: Both `times` and `new_times` must be non-negative, unique, INCREASING numbers
// Assumptions (we always use the end time point in an interval):
// 1) x(t < t1) = x(t1)
// 2) x(t_k < t < t_{k+1}) = x(t_{k+1})
// 3) x(t > t_max) = x(t_max)
// [[Rcpp::export]]
NumericMatrix rcpp_wcumsum_mat_interp(const NumericMatrix& x,
                                      const NumericVector& times,
                                      const NumericVector& new_times) {
  int n_rows = x.nrow();
  int n_times = x.ncol(); // Number of original time points
  int n_times_new = new_times.length(); // Number of new time points
  NumericMatrix mat(n_rows, n_times_new);

  // Compute time differences (weights)
  std::vector<double> delta_t(n_times);
  delta_t[0] = times[0]; // First interval from time 0
  for (int j = 1; j < n_times; j++) {
    delta_t[j] = times[j] - times[j - 1];
  }

  // Process each row
  for (int i = 0; i < n_rows; i++) {
    double cum_sum = 0.0;  // Cumulative sum for the row
    int j = 0;             // Pointer for `times`

    for (int k = 0; k < n_times_new; k++) {
      double t_new = new_times[k];

      // If `t_new` is before the first time point, take the rectangle area only
      if (t_new < times[0]) {
        mat(i, k) = x(i, 0) * t_new;
        continue;
      }

      // Accumulate full intervals up to `t_new`
      while (j < n_times && times[j] <= t_new) {
        cum_sum += x(i, j) * delta_t[j];
        j++;
      }

      // Add missing rectangle area
      double extra_area = (j == n_times)
        // `t_new` > tmax
        ? x(i, n_times - 1) * (t_new - times[n_times - 1])
        // `t_new` in between, j points to the next time
        : x(i, j) * (t_new - times[j - 1]);

      mat(i, k) = cum_sum + extra_area;
    }
  }

  return mat;
}
