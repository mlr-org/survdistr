#include <Rcpp.h>
using namespace Rcpp;

// Per-row interpolation of a matrix with survival values
// * `x`: survival matrix [observations x times]
// * `times`: original time points (increasing, unique, non-negative)
// * `eval_times`: new time points to interpolate (increasing, unique, non-negative)
// * `constant`: if true, stepwise constant (left-continuous); otherwise linear interpolation
// * type: "surv", "cdf", or "cif"
// [[Rcpp::export]]
NumericMatrix c_mat_interp(const NumericMatrix& x,
                           const NumericVector& times,
                           const NumericVector& eval_times,
                           bool constant = true,
                           const std::string& type = "surv") {
  // x => [obs x times]
  int n_rows = x.nrow(); // observations
  int n_times = x.ncol(); // original time points
  int n_times_eval = eval_times.length(); // requested time points

  // Baseline value for extrapolation
  double BASE_DEFAULT = (type == "surv") ? 1.0 : 0.0;

  NumericMatrix mat(n_rows, n_times_eval);

  // Iterate over observations
  for (int i = 0; i < n_rows; i++) {
    int j = 0; // index for scanning times

    // Iterate over requested time points
    for (int k = 0; k < n_times_eval; k++) {
      double t_new = eval_times[k];

      // extrapolation before first time
      if (t_new < times[0]) {
        if (constant) {
          mat(i, k) = BASE_DEFAULT;
        } else {
          if (times[0] == 0.0) {
            mat(i, k) = x(i, 0); // avoid div-by-zero
          } else {
            mat(i, k) = BASE_DEFAULT + t_new * (x(i, 0) - BASE_DEFAULT) / times[0];
          }
        }
        continue;
      }

      // extrapolation after last time
      if (t_new >= times[n_times - 1]) {
        if (constant || n_times <= 1) {
          // Use last value for constant extrapolation (left-continuity)
          mat(i, k) = x(i, n_times - 1);
        } else {
          // Linear extrapolation using the last time interval
          double t1 = times[n_times - 2], t2 = times[n_times - 1];
          double x1 = x(i, n_times - 2), x2 = x(i, n_times - 1);
          double val = x2 + (t_new - t2) * (x2 - x1) / (t2 - t1);

          // Adjust bounds
          if (type == "surv") {
            // avoid S(t) < 0
            mat(i, k) = std::max(0.0, val);
          } else {
            // avoid CDF(t) > 1 or CIF(t) > 1
            mat(i, k) = std::min(1.0, val);
          }
        }
        continue;
      }

      // advance j until t_new is in [times[j], times[j+1])
      while (j < n_times - 1 && t_new >= times[j + 1]) j++;

      // t_new in [times[j], times[j+1])
      if (constant) {
        // Constant interpolation (use value at the lower bound)
        mat(i, k) = x(i, j);
      } else {
        // Linear interpolation
        double t1 = times[j], t2 = times[j + 1];
        double x1 = x(i, j), x2 = x(i, j + 1);
        mat(i, k) = x1 + (t_new - t1) * (x2 - x1) / (t2 - t1);
      }
    }
  }

  return mat;
}

// Interpolation of a single survival/CDF/CIF curve/vector
// [[Rcpp::export]]
NumericVector c_vec_interp(const NumericVector& x,
                           const NumericVector& times,
                           const NumericVector& eval_times,
                           bool constant = true,
                           const std::string& type = "surv") {
  int n_times = x.size();
  int n_times_eval = eval_times.size();
  NumericVector out(n_times_eval);

  double BASE_DEFAULT = (type == "surv") ? 1.0 : 0.0;

  int j = 0; // index in original times

  for (int k = 0; k < n_times_eval; ++k) {
    double t_new = eval_times[k];

    // extrapolation before first time
    if (t_new < times[0]) {
      if (constant) {
        out[k] = BASE_DEFAULT;
      } else {
        if (times[0] == 0.0) {
          out[k] = x[0];
        } else {
          out[k] = BASE_DEFAULT + t_new * (x[0] - BASE_DEFAULT) / times[0];
        }
      }
      continue;
    }

    // extrapolation after last time
    if (t_new >= times[n_times - 1]) {
      if (constant || n_times <= 1) {
        out[k] = x[n_times - 1];
      } else {
        double t1 = times[n_times - 2], t2 = times[n_times - 1];
        double x1 = x[n_times - 2], x2 = x[n_times - 1];
        double val = x2 + (t_new - t2) * (x2 - x1) / (t2 - t1);

        // Adjust bounds
        if (type == "surv") {
          out[k] = std::max(0.0, val);
        } else {
          out[k] = std::min(1.0, val);
        }
      }
      continue;
    }

    // advance j until t_new is in [times[j], times[j+1])
    while (j < n_times - 1 && t_new >= times[j + 1]) j++;

    if (constant) {
      out[k] = x[j];
    } else {
      double t1 = times[j], t2 = times[j + 1];
      double x1 = x[j], x2 = x[j + 1];
      out[k] = x1 + (t_new - t1) * (x2 - x1) / (t2 - t1);
    }
  }

  return out;
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
