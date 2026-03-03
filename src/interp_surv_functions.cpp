#include <Rcpp.h>
using namespace Rcpp;

enum InterpMethod {
  CONST_SURV,
  LINEAR_SURV,
  CONST_HAZ
};

InterpMethod parse_method(const std::string& method) {
  if (method == "const_surv")  return CONST_SURV;
  if (method == "linear_surv") return LINEAR_SURV;
  if (method == "const_haz")   return CONST_HAZ;
  stop("Unknown interpolation method.");
}

// [[Rcpp::export]]
NumericMatrix c_interp_surv_mat(
  const NumericMatrix& x, // survival matrix [obs x times]
  const NumericVector& times, // anchor time points (increasing, unique, non-negative)
  const NumericVector& eval_times, // new time points to interpolate (increasing, unique, non-negative)
  const std::string& method = "const_surv" // interpolation method
) {
  const int n_rows = x.nrow(); // observations
  const int B = times.size(); // number of anchor points (t_1,...,t_B)
  const int n_eval = eval_times.size(); // requested time points

  InterpMethod m = parse_method(method);
  NumericMatrix result(n_rows, n_eval);

  // Iterate over observations
  for (int i = 0; i < n_rows; i++) {
    // index of largest anchor <= t (j == -1 means t < times[0] => before first anchor)
    int j = -1;

    // Iterate over requested time points
    for (int k = 0; k < n_eval; k++) {
      double t = eval_times[k];

      // t = 0 => S(0) = 1
      if (t == 0.0) {
        result(i, k) = 1.0;
        continue;
      }

      // Advance j so times[j] <= t < times[j+1] (if possible)
      while (j < B-1 && times[j + 1] <= t) {
        ++j;
      }

      // Exact anchor
      if (j >= 0 && times[j] == t) {
        result(i, k) = x(i, j);
        continue;
      }

      // after last anchor: t > times[B-1] (extrapolation)
      if (j == B - 1) {
        if (m == CONST_SURV) {
          result(i, k) = x(i, B - 1);
        } else if (m == LINEAR_SURV) {
          // choose last interval; if only one anchor use (0, times[0])
          double t_left, t_right, S_left, S_right;

          if (B == 1) {
            t_left = 0.0;
            t_right = times[0];
            S_left = 1.0;
            S_right = x(i, 0);
          } else {
            t_left = times[B - 2];
            t_right = times[B - 1];
            S_left = x(i, B - 2);
            S_right = x(i, B - 1);
          }
          double slope = (S_right - S_left) / (t_right - t_left);
          double val = S_right + slope * (t - t_right);
          result(i, k) = std::max(0.0, std::min(1.0, val));
        } else { // CONST_HAZ placeholder
          result(i, k) = x(i, B - 1);
        }
        continue;
      }

      // t in (t_j, t_{j+1}) where j = -1 means we are at (0, times[0])
      double t_left, t_right, S_left, S_right;
      if (j == -1) {
        t_left = 0.0;
        t_right = times[0];
        S_left = 1.0;
        S_right = x(i, 0);
      } else {
        t_left = times[j];
        t_right = times[j + 1];
        S_left = x(i, j);
        S_right = x(i, j + 1);
      }

      if (m == CONST_SURV) {
        result(i, k) = S_left;
      } else if (m == LINEAR_SURV) {
        double val = S_left + (t - t_left) * (S_right - S_left) / (t_right - t_left);
        result(i, k) = std::max(0.0, std::min(1.0, val));
      } else { // CONST_HAZ placeholder
        result(i, k) = S_left;
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
NumericVector c_interp_surv_vec(
  const NumericVector& x, // survival vector [times]
  const NumericVector& times, // anchor times (increasing, unique, non-negative)
  const NumericVector& eval_times, // new time points to interpolate (increasing, unique, non-negative)
  const std::string& method = "const_surv" // interpolation method
) {
  const int B = times.size(); // number of anchor points (t_1,...,t_B)
  const int n_eval = eval_times.size(); // requested time points

  InterpMethod m = parse_method(method);
  NumericVector result(n_eval);

  // index of largest anchor <= t (j == -1 means t < times[0] => before first anchor)
  int j = -1;

  // Iterate over requested time points
  for (int k = 0; k < n_eval; k++) {
    double t = eval_times[k];

    // t = 0 => S(0) = 1
    if (t == 0.0) {
      result[k] = 1.0;
      continue;
    }

    // Advance j so times[j] <= t < times[j+1] (if possible)
    while (j < B-1 && times[j + 1] <= t) {
      ++j;
    }

    // Exact anchor
    if (j >= 0 && times[j] == t) {
      result[k] = x[j];
      continue;
    }

    // after last anchor: t > times[B-1] (extrapolation)
    if (j == B - 1) {
      if (m == CONST_SURV) {
        result[k] = x[B - 1];
      } else if (m == LINEAR_SURV) {
        // choose last interval; if only one anchor use (0, times[0])
        double t_left, t_right, S_left, S_right;

        if (B == 1) {
          t_left = 0.0;
          t_right = times[0];
          S_left = 1.0;
          S_right = x[0];
        } else {
          t_left = times[B - 2];
          t_right = times[B - 1];
          S_left = x[B - 2];
          S_right = x[B - 1];
        }
        double slope = (S_right - S_left) / (t_right - t_left);
        double val = S_right + slope * (t - t_right);
        result[k] = std::max(0.0, std::min(1.0, val));
      } else { // CONST_HAZ placeholder
        result[k] = x[B - 1];
      }
      continue;
    }

    // t in (t_j, t_{j+1}) where j = -1 means we are at (0, times[0])
    double t_left, t_right, S_left, S_right;
    if (j == -1) {
      t_left = 0.0;
      t_right = times[0];
      S_left = 1.0;
      S_right = x[0];
    } else {
      t_left = times[j];
      t_right = times[j + 1];
      S_left = x[j];
      S_right = x[j + 1];
    }

    if (m == CONST_SURV) {
      result[k] = S_left;
    } else if (m == LINEAR_SURV) {
      double val = S_left + (t - t_left) * (S_right - S_left) / (t_right - t_left);
      result[k] = std::max(0.0, std::min(1.0, val));
    } else { // CONST_HAZ placeholder
      result[k] = S_left;
    }
  }
  
  return result;
}

// [[Rcpp::export]]
NumericMatrix c_interp_cif_mat(
  const NumericMatrix& x, // CIF matrix [obs x times]
  const NumericVector& times, // anchor time points (increasing, unique, non-negative)
  const NumericVector& eval_times // new time points to interpolate (increasing, unique, non-negative)
) {
  int n_rows = x.nrow();
  int n_times = x.ncol();
  int n_times_eval = eval_times.length();

  const double BASE_DEFAULT = 0.0;

  NumericMatrix mat(n_rows, n_times_eval);

  for (int i = 0; i < n_rows; i++) {

    int j = 0;

    for (int k = 0; k < n_times_eval; k++) {

      double t_new = eval_times[k];

      if (t_new < times[0]) {
        mat(i, k) = BASE_DEFAULT;
        continue;
      }

      if (t_new >= times[n_times - 1]) {
        mat(i, k) = x(i, n_times - 1);
        continue;
      }

      while (j < n_times - 1 && t_new >= times[j + 1])
        j++;

      mat(i, k) = x(i, j);
    }
  }

  return mat;
}

// [[Rcpp::export]]
NumericVector c_mat_interp_pointwise(const NumericMatrix& x,
                                     const NumericVector& times,
                                     const NumericVector& eval_times,
                                     bool constant = true,
                                     const std::string& type = "surv") {
  int n_obs = x.nrow();
  int n_times = x.ncol();

  if (eval_times.size() != n_obs) {
    stop("Length of 'eval_times' must match number of rows of 'x'.");
  }

  NumericVector out(n_obs);
  double BASE_DEFAULT = (type == "surv") ? 1.0 : 0.0;

  for (int i = 0; i < n_obs; ++i) {
    double t_new = eval_times[i];

    // extrapolation before first time
    if (t_new < times[0]) {
      if (constant) {
        out[i] = BASE_DEFAULT;
      } else {
        if (times[0] == 0.0) {
          out[i] = x(i, 0);
        } else {
          out[i] = BASE_DEFAULT + t_new * (x(i, 0) - BASE_DEFAULT) / times[0];
        }
      }
      continue;
    }

    // extrapolation after last time
    if (t_new >= times[n_times - 1]) {
      if (constant || n_times <= 1) {
        // Use last value for constant extrapolation (left-continuity)
        out[i] = x(i, n_times - 1);
      } else {
        // Linear extrapolation using the last time interval
        double t1 = times[n_times - 2], t2 = times[n_times - 1];
        double x1 = x(i, n_times - 2), x2 = x(i, n_times - 1);
        double val = x2 + (t_new - t2) * (x2 - x1) / (t2 - t1);
        out[i] = (type == "surv") ? std::max(0.0, val) : std::min(1.0, val);
      }
      continue;
    }

    // interpolation
    // advance j until t_new is in [times[j], times[j+1])
    int j = 0;
    while (j < n_times - 1 && t_new >= times[j + 1]) j++;

    if (constant) {
      // Constant interpolation (use value at the lower bound)
      out[i] = x(i, j);
    } else {
      // Linear interpolation
      double t1 = times[j], t2 = times[j + 1];
      double x1 = x(i, j), x2 = x(i, j + 1);
      out[i] = x1 + (t_new - t1) * (x2 - x1) / (t2 - t1);
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
