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
  const int n_obs = x.nrow(); // observations
  const int B = times.size(); // number of anchor points (t_1,...,t_B) - (must equal x.ncol())
  const int n_eval = eval_times.size(); // requested time points

  InterpMethod m = parse_method(method);
  NumericMatrix result(n_obs, n_eval);

  // Iterate over observations
  for (int i = 0; i < n_obs; i++) {
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

      // Advance j to the interval containing t: times[j] ≤ t < times[j+1]
      while (j < B - 1 && t >= times[j + 1]) {
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
