#include <Rcpp.h>
using namespace Rcpp;

// for mapping h(t) or f(t), discrete or continuous output, to S(t),
// depending on the interpolation method

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
