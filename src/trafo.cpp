#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
// Map row-wise h(t) or f(t), discrete or continuous, to S(t) at given anchor time points

// Discrete density => survival: S_j = 1 - sum_{k=1}^j f_k
// [[Rcpp::export]]
NumericMatrix c_disc_dens_to_surv_mat(const NumericMatrix& x) {
  const int n_rows = x.nrow();
  const int n_times = x.ncol();
  NumericMatrix surv(n_rows, n_times);

  for (int i = 0; i < n_rows; i++) {
    double cum_sum = 0.0;
    for (int j = 0; j < n_times; j++) {
      cum_sum += x(i, j);
      surv(i, j) = 1.0 - cum_sum;
    }
  }

  return surv;
}

// Discrete hazard => survival: S_j = prod_{k=1}^j (1 - h_k)
// [[Rcpp::export]]
NumericMatrix c_disc_haz_to_surv_mat(const NumericMatrix& x) {
  const int n_rows = x.nrow();
  const int n_times = x.ncol();
  NumericMatrix surv(n_rows, n_times);

  for (int i = 0; i < n_rows; i++) {
    double prod_surv = 1.0;
    for (int j = 0; j < n_times; j++) {
      prod_surv *= (1.0 - x(i, j));
      surv(i, j) = prod_surv;
    }
  }

  return surv;
}

// Continuous density => survival: S_j = 1 - sum_{k=1}^j f_k * Delta_k
// [[Rcpp::export]]
NumericMatrix c_cont_dens_to_surv_mat(
  const NumericMatrix& x,
  const NumericVector& times
) {
  const int n_rows = x.nrow();
  const int n_times = x.ncol();
  NumericMatrix surv(n_rows, n_times);

  std::vector<double> delta_t(n_times);
  delta_t[0] = times[0];
  for (int j = 1; j < n_times; j++) {
    delta_t[j] = times[j] - times[j - 1];
  }

  for (int i = 0; i < n_rows; i++) {
    double cum_int = 0.0;
    for (int j = 0; j < n_times; j++) {
      cum_int += x(i, j) * delta_t[j];
      surv(i, j) = 1.0 - cum_int;
    }
  }

  return surv;
}

// Continuous hazards => survival: S_j = exp(-sum_{k=1}^j lambda_k * Delta_k)
// [[Rcpp::export]]
NumericMatrix c_cont_haz_to_surv_mat(
  const NumericMatrix& x,
  const NumericVector& times
) {
  const int n_rows = x.nrow();
  const int n_times = x.ncol();
  NumericMatrix surv(n_rows, n_times);

  std::vector<double> delta_t(n_times);
  delta_t[0] = times[0];
  for (int j = 1; j < n_times; j++) {
    delta_t[j] = times[j] - times[j - 1];
  }

  for (int i = 0; i < n_rows; i++) {
    double cum_haz = 0.0;
    for (int j = 0; j < n_times; j++) {
      cum_haz += x(i, j) * delta_t[j];
      surv(i, j) = std::exp(-cum_haz);
    }
  }

  return surv;
}
