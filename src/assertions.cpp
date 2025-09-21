#include <Rcpp.h>
using namespace Rcpp;

// Per-row validation of survival / CDF / CIF matrices
// `x`: non-empty, numeric matrix [obs x times]
// `type`: "surv", "cdf", or "cif"
// [[Rcpp::export]]
bool c_assert_prob_matrix(const NumericMatrix& x, const std::string& type = "surv") {
  int nrows = x.nrow();
  int ncols = x.ncol();

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      double value = x(i, j);

      // Probabilities must be in [0, 1]
      if (value < 0.0 || value > 1.0) {
        return false;
      }

      if (type == "surv") {
        // survival should be non-increasing
        if (j > 0 && value > x(i, j - 1)) return false;
      } else {
        // CDF / CIF should be non-decreasing
        if (j > 0 && value < x(i, j - 1)) return false;
      }
    }
  }

  return true;
}
