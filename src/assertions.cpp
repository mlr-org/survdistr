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

// Clamp survival values into [eps, 1] with warning if large violations
// [[Rcpp::export]]
NumericMatrix c_clamp_surv(
  NumericMatrix surv,
  const double eps = 1e-12,
  const double tol = 1e-10
) {
  const int n = surv.size();
  bool warn_flag = false;

  for (int k = 0; k < n; ++k) {
    double val = surv[k];

    // detect meaningful violation
    if (val < -tol || val > 1.0 + tol) {
      warn_flag = true;
    }

    // clamp
    if (val > 1.0) {
      surv[k] = 1.0;
    } else if (val < eps) {
      surv[k] = eps;
    }
  }

  if (warn_flag) {
    Rcpp::warning(
      "Survival values outside [0,1] beyond tolerance; values were clamped."
    );
  }

  return surv;
}
