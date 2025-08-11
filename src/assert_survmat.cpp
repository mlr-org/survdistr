#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
bool rcpp_assert_surv_matrix(const NumericMatrix& mat) {
  int nrows = mat.nrow();
  int ncols = mat.ncol();

  // Handle edge cases for empty matrices
  if (nrows == 0 || ncols == 0) {
    stop("Input matrix must have at least one row and one column.");
  }

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      double value = mat(i, j);

      // S(t) in [0, 1]
      if (value < 0.0 || value > 1.0) {
        return false;
      }

      // S(t) > S(t-1) is not acceptable by definition
      if (j > 0 && value > mat(i, j - 1)) {
        return false;
      }
    }
  }

  return true;
}
