#' @title Assert survival matrix
#'
#' @description Asserts if the given input is a (discrete) survival
#' probabilities matrix using [Rcpp] code.
#' The following checks are performed:
#'
#' 1. All values are probabilities, i.e. \eqn{S(t) \in [0,1]}
#' 2. Column names correspond to time points and should therefore be coercable to
#' `numeric` and increasing
#' 3. Per row (observation), the survival probabilities decrease non-strictly, i.e.
#' \eqn{S(t) \ge S(t+1)}
#'
#' @param x (`matrix()`)\cr
#' A survival probability matrix.
#' Rows are observations and columns are (increasing) time points.
#'
#' @return if the assertion fails an error occurs, otherwise `NULL` is returned
#' invisibly.
#'
#' @examples
#' x = matrix(data = c(1,0.6,0.4,0.8,0.8,0.7), nrow = 2, ncol = 3, byrow = TRUE)
#' colnames(x) = c(12, 34, 42)
#' x
#'
#' assert_surv_matrix(x)
#'
#' @export
assert_surv_matrix = function(x) {
  assert_matrix(x, any.missing = FALSE, min.rows = 1, min.cols = 1, col.names = "named")
  times = assert_numeric(as.numeric(colnames(x)), lower = 0, unique = TRUE, sorted = TRUE,
                         any.missing = FALSE, null.ok = FALSE)

  if (!rcpp_assert_surv_matrix(x)) {
    stop("Survival probabilities must be (non-strictly) decreasing and between [0, 1]")
  }

  if (times[1] == 0 && !all(x[, 1] == 1)) {
    stop("First time point is equal to zero and S(t = 0) != 1")
  }

  invisible(NULL)
}

#' @title Assert hazard matrix
#'
#' @description Asserts if the given input is a (discrete) hazard matrix.
#' The following checks are performed:
#'
#' 1. All values are non-negative, i.e. \eqn{h(t) \ge 0}
#' 2. Column names correspond to time points and should therefore be coercable to
#' `numeric` and increasing
#'
#' @param x (`matrix()`)\cr
#' A matrix of hazards.
#' Rows are observations and columns are (increasing) time points.
#'
#' @return if the assertion fails an error occurs, otherwise `NULL` is returned
#' invisibly.
#'
#' @examples
#' x = matrix(data = c(0,0.2,0.1,0.5,0.05,0.7), nrow = 2, ncol = 3, byrow = TRUE)
#' colnames(x) = c(12, 34, 42)
#' x
#'
#' assert_hazard_matrix(x)
#'
#' @export
assert_hazard_matrix = function(x) {
  assert_matrix(x, any.missing = FALSE, min.rows = 1, min.cols = 1, col.names = "named")
  assert_numeric(as.numeric(colnames(x)), unique = TRUE, sorted = TRUE,
                 any.missing = FALSE, null.ok = FALSE)
  stopifnot(x >= 0)

  invisible(NULL)
}
