#include <Rcpp.h>
#include "interp_helpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix c_interp_density_mat(
  const NumericMatrix& x, // survival matrix [obs x times]
  const NumericVector& times, // anchor time points (increasing, unique, non-negative)
  const NumericVector& eval_times, // new time points to interpolate (increasing, unique, non-negative)
  const std::string& method = "const_surv" // interpolation method
) {
  const int n_obs = x.nrow(); // observations
  const int B = times.size(); // number of anchor points (t_1,...,t_B) - (must equal x.ncol())
  const int n_eval = eval_times.size(); // requested time points
  const double eps = 1e-12; // lower value for zero survival values under constant hazard interpolation

  InterpMethod m = parse_method(method);
  NumericMatrix density(n_obs, n_eval);

  // Iterate over observations
  for (int i = 0; i < n_obs; i++) {
    // index of largest anchor <= t (j == -1 means t < times[0] => before first anchor)
    int j = -1;

    // Iterate over requested time points
    for (int k = 0; k < n_eval; k++) {
      double t = eval_times[k];

      // ----- t = 0 -----
      // t = 0 => f(0) = 0 if only one anchor at 0, else the value from the first interval (0, times[0])
      if (t == 0.0) {
        // density is zero except at anchors >0
        if (m == CONST_SURV) {
          density(i, k) = 0.0;
          continue;
        }

        // Pre-initialize t_left and S_left
        double t_left = 0.0, S_left = 1.0, t_right, S_right;

        if (times[0] > 0.0) {
          // first interval is (0, t1)
          t_right = times[0];
          S_right = x(i,0);
        } else if (B > 1) {
          // first interval is (t1 = 0, t2), i.e first anchor is at 0
          t_right = times[1];
          S_right = x(i,1);
        } else {
          // only anchor at t = 0 => no interval information
          density(i, k) = 0.0;
          continue;
        }

        double delta = t_right - t_left;

        if (m == CONST_DENS) {
          density(i, k) = (S_left - S_right) / delta;
        } else { // CONST_HAZ
          density(i, k) = -std::log(S_right > 0.0 ? S_right : eps) / delta;
        }
        continue;
      }

      // Advance j to the interval containing t: times[j] ≤ t < times[j+1]
      while (j < B-1 && t >= times[j+1]) {
        ++j;
      }

      // ----- Exact anchor (t > 0) -----
      if (j >= 0 && times[j] == t) {
        double t_left, t_right, S_left, S_right;

        if (j == 0) {
          t_left = 0.0;
          t_right = times[0];
          S_left = 1.0;
          S_right = x(i, 0);
        } else {
          t_left = times[j-1];
          t_right = times[j];
          S_left = x(i, j-1);
          S_right = x(i, j);
        }
        double delta = t_right - t_left;

        if (m == CONST_SURV) {
          // f(t) != 0 only at anchors for constant survival interpolation
          density(i, k) = S_left - S_right;
        } else if (m == CONST_DENS) {
          density(i, k) = (S_left - S_right) / delta;
        } else { // CONST_HAZ
          // avoid unnessary division by zero
          if (S_left == 0.0) {
            density(i, k) = 0.0;
          } else {
            double ratio = S_right <= 0.0 ? std::min(eps, S_left) / S_left : S_right / S_left;
            double haz = -std::log(ratio) / delta;
            density(i, k) = haz * S_right;
          }
        }
        continue;
      }

      // ----- Extrapolation => after last anchor: t > times[B-1] -----
      if (j == B - 1) {
        // handle constant survival extrapolation early on for speed
        if (m == CONST_SURV) {
          density(i, k) = 0.0;
          continue;
        }

        // choose last interval; if only one anchor use (0, times[0])
        double t_left, t_right, S_left, S_right;

        if (B == 1) {
          t_left = 0.0;
          t_right = times[0];
          S_left = 1.0;
          S_right = x(i, 0);
        } else {
          t_left = times[B-2];
          t_right = times[B-1];
          S_left = x(i, B-2);
          S_right = x(i, B-1);
        }

        // no need to extrapolate if survival is already 0 at last anchor
        if (S_right <= 0.0) {
          density(i, k) = 0.0;
          continue;
        }

        double delta = t_right - t_left;

        // degenerate interval (e.g. single anchor at t = 0)
        if (delta == 0.0) {
          density(i, k) = 0.0;
          continue;
        }

        if (m == CONST_DENS) {
          double fB = (S_left - S_right) / delta;
          double t_star = t_right + S_right / fB;

          // If S_left == S_right > 0 (last interval is flat) then fB = 0 and t_star = +Inf => f(t) = fB = 0
          density(i,k) = (t < t_star) ? fB : 0.0;
        } else { // CONST_HAZ
          // definitely here S_left > 0 and S_right > 0 since we handle the S_right == 0 case above
          double alpha = (t - t_right) / delta;
          double ratio = S_right / S_left;
          double surv = S_right * std::pow(ratio, alpha);
          double haz = -std::log(ratio) / delta;
          density(i, k) = haz * surv;
        }
        continue;
      }

      // ----- Interpolation => t in (t_j, t_{j+1}) where j = -1 means we are at (0, times[0] > 0) -----
      // handle constant survival interpolation early on for speed
      if (m == CONST_SURV) {
        density(i, k) = 0.0;
        continue;
      }

      double t_left, t_right, S_left, S_right;
      if (j == -1) {
        t_left = 0.0;
        t_right = times[0];
        S_left = 1.0;
        S_right = x(i, 0);
      } else {
        t_left = times[j];
        t_right = times[j+1];
        S_left = x(i, j);
        S_right = x(i, j+1);
      }
      double delta = t_right - t_left;

      if (m == CONST_DENS) {
        density(i, k) = (S_left - S_right) / delta;
      } else { // CONST_HAZ
        // avoid unnessary division by zero
        if (S_left == 0.0) {
          density(i, k) = 0.0;
        } else {
          double alpha = (t - t_left) / delta;
          double ratio = S_right <= 0.0 ? std::min(eps, S_left) / S_left : S_right / S_left;
          double surv = S_left * std::pow(ratio, alpha);
          double haz = -std::log(ratio) / delta;
          density(i, k) = haz * surv;
        }
      }
    }
  }

  return density;
}
