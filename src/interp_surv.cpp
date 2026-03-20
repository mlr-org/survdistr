#include <Rcpp.h>
#include "interp_helpers.h"
using namespace Rcpp;

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
  NumericMatrix survival(n_obs, n_eval);

  // Iterate over observations
  for (int i = 0; i < n_obs; i++) {
    // index of largest anchor <= t (j == -1 means t < times[0] => before first anchor)
    int j = -1;

    // Iterate over requested time points
    for (int k = 0; k < n_eval; k++) {
      double t = eval_times[k];

      // ----- t = 0 -----
      // t = 0 => S(0) = 1
      if (t == 0.0) {
        survival(i, k) = 1.0;
        continue;
      }

      // Advance j to the interval containing t: times[j] ≤ t < times[j+1]
      while (j < B - 1 && t >= times[j + 1]) {
        ++j;
      }

      // ----- Exact anchor (t > 0) -----
      if (j >= 0 && times[j] == t) {
        survival(i, k) = x(i, j);
        continue;
      }

      // ----- Extrapolation => after last anchor: t > times[B-1] -----
      if (j == B - 1) {
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

        // no need to extrapolate if survival is already 0 at last anchor
        if (S_right <= 0.0) {
          survival(i, k) = 0.0;
          continue;
        }

        // handle constant survival extrapolation early on for speed
        if (m == CONST_SURV) {
          survival(i, k) = S_right;
          continue;
        }

        double delta = t_right - t_left;

        // degenerate interval (e.g. single anchor at t = 0)
        if (delta == 0.0) {
          survival(i, k) = S_right;
          continue;
        }

        if (m == CONST_DENS) {
          double fB = (S_left - S_right) / delta;
          double t_star = t_right + S_right / fB;

          // If S_left == S_right > 0 (last interval is flat) then fB = 0 and t_star = +Inf => S(t) = S_right
          if (t >= t_star) {
            survival(i, k) = 0.0;
          } else {
            survival(i, k) = S_right - fB * (t - t_right);
          }
        } else { // CONST_HAZ
          // definitely here S_left > 0 since it's been handled above with the S_right = 0 case
          double alpha = (t - t_right) / delta;
          double ratio = S_right / S_left;
          double val = S_right * std::pow(ratio, alpha);

          survival(i, k) = std::max(0.0, std::min(1.0, val));
        }
        continue;
      }

      // ----- Interpolation => t in (t_j, t_{j+1}) where j = -1 means we are at (0, times[0] > 0) -----
      // handle constant survival interpolation early on for speed
      double S_left = (j == -1) ? 1.0 : x(i, j);
      if (m == CONST_SURV) {
        survival(i, k) = S_left;
        continue;
      }

      double t_left, t_right, S_right;
      if (j == -1) {
        t_left = 0.0;
        t_right = times[0];
        S_right = x(i, 0);
      } else {
        t_left = times[j];
        t_right = times[j + 1];
        S_right = x(i, j + 1);
      }
      double delta = t_right - t_left;
      double alpha = (t - t_left) / delta;

      if (m == CONST_DENS) {
        double val = S_left + alpha * (S_right - S_left);
        survival(i, k) = std::max(0.0, std::min(1.0, val));
      } else { // CONST_HAZ
        // avoid unnecessary division by zero
        if (S_left == 0.0) {
          survival(i, k) = 0.0;
        } else {
          double ratio = S_right / S_left;
          double val = S_left * std::pow(ratio, alpha);
          survival(i, k) = std::max(0.0, std::min(1.0, val));
        }
      }
    }
  }

  return survival;
}
