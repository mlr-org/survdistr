#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix c_interp_cif_mat(
  const NumericMatrix& x, // CIF matrix [obs x times]
  const NumericVector& times, // anchor time points (increasing, unique, non-negative)
  const NumericVector& eval_times // new time points to interpolate (increasing, unique, non-negative)
) {
  const int n_obs = x.nrow();
  const int B = times.size(); // number of anchor points (must equal x.ncol())
  const int n_eval = eval_times.size();

  NumericMatrix result(n_obs, n_eval);

  for (int i = 0; i < n_obs; i++) {
    int j = 0; // index of largest anchor <= current eval_time

    for (int k = 0; k < n_eval; k++) {
      double t = eval_times[k];

      // Before first anchor: CIF = 0
      if (t < times[0]) {
        result(i, k) = 0.0;
        continue;
      }

      // After last anchor: constant extrapolation (use last value)
      if (t >= times[B - 1]) {
        result(i, k) = x(i, B - 1);
        continue;
      }

      // Advance j to the interval containing t: times[j] ≤ t < times[j+1]
      while (j < B - 1 && t >= times[j + 1]) {
        ++j;
      }

      // Constant interpolation: value at left anchor
      result(i, k) = x(i, j);
    }
  }

  return result;
}
