#include <Rcpp.h>
using namespace Rcpp;

// Validation (rowwise) of survival, hazard, density, CDF and CIF probaiblity matrices
// [[Rcpp::export]]
bool c_assert_prob_matrix(
  const NumericMatrix& x,
  const std::string& type = "surv"
) {
  const int nrows = x.nrow();
  const int ncols = x.ncol();

  bool check_inc = false;
  bool check_dec = false;

  if (type == "surv") {
    check_dec = true;
  } else if (type == "cdf" || type == "cif") {
    check_inc = true;
  } else if (type == "prob" || type == "haz" || type == "dens") {
    // only range check, i.e. probabilities in [0, 1]
  } else {
    stop("Unknown type.");
  }

  for (int i = 0; i < nrows; ++i) {
    double prev = x(i, 0);
    if (prev < 0.0 || prev > 1.0) return false;

    for (int j = 1; j < ncols; ++j) {
      double val = x(i, j);

      if (val < 0.0 || val > 1.0) return false;
      if (check_dec && val > prev) return false;
      if (check_inc && val < prev) return false;
      prev = val;
    }
  }

  return true;
}
