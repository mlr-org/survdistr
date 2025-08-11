#' @title Constructor for survDistr
#'
#' @description
#' Creates a [survDistr] object, a specialized container designed for storing
#' and managing prediction outputs from survival models in single-event settings.
#' This includes models such as Cox proportional hazards, random survival forests,
#' and other classical or machine learning-based survival estimators.
#'
#' The main prediction data type can be a survival or a hazard matrix, where
#' **rows represent observations and columns represent time points**.
#'
#' @details
#' The input matrix (survival probabilities \eqn{S(t)} or hazard \eqn{h(t)})
#' is stored in the `$data` slot, while the interpolation type for subsequent S3
#' methods is stored in the `$inter_type` slot.
#'
#' During construction, the functions [assert_surv_matrix()] or [assert_hazard_matrix()]
#' are used to validate the input data matrix according to the given `data_type`.
#'
#' @param x (`matrix()`)\cr
#' A numeric matrix of either survival probabilities (values between 0 and 1) or
#' hazard values (non-negative values).
#' Column names must correspond to time points.
#' @template param_data_type
#' @template param_inter_type
#' @param ... currently not used
#'
#' @section Available S3 Methods:
#' - [print()]
#' - [times()]
#' - [survival()]
#' - [cumhazard()]
#' - [hazard()]
#' - [cdf()]
#' - [pdf()]
#'
#' @return An object of class [survDistr] containing the survival matrix.
#' @examples
#' x = matrix(data = c(1,0.6,0.4,0.8,0.8,0.7), nrow = 2,
#'            ncol = 3, byrow = TRUE)
#' colnames(x) = c(12, 34, 42)
#' x = survDistr(x)
#' x
#'
#' # stored survival matrix
#' x$data
#'
#' # type of interpolation to use in S3 methods
#' x$inter_type
#'
#' # time points
#' times(x)
#'
#' # S(t) at given time points (constant interpolation)
#' survival(x, times = c(10, 30, 42, 45))
#' # same but with linear interpolation
#' x$inter_type = "linear_surv"
#' survival(x, times = c(10, 30, 42, 45))
#' # time points can be unordered and duplicated
#' survival(x, times = c(10,30,10,50))
#'
#' # Cumulative hazard
#' cumhazard(x)
#'
#' # hazard
#' hazard(x)
#'
#' # cumulative distribution function
#' cdf(x)
#'
#' # probability density function
#' pdf(x)
#'
#' @export
survDistr = function(x, data_type = "surv", inter_type = "const_surv", ...) {
  assert_choice(data_type, c("surv", "haz"))
  if (data_type == "surv") {
    assert_surv_matrix(x)
  } else {
    assert_hazard_matrix(x)
  }

  assert_choice(inter_type, c("const_surv", "linear_surv", "const_haz"))

  # checks for data_type and inter_type
  if (data_type == "haz" && inter_type == "linear_surv")  {
    stop("Hazard metric and piece-wise linear survival interpolation is not supported.")
  }

  # Create the S3 object
  structure(
    list(data = x, data_type = data_type, inter_type = inter_type),
    class = "survDistr"
  )
}

#' @title Print method for survDistr
#' @description
#' Displays summary information about a [survDistr] object, including
#' the number of observations and time points.
#' @template param_x
#' @export
print.survDistr = function(x, ...) {
  nrows = nrow(x$data)
  ncols = ncol(x$data)
  cat("A [", nrows, " x ", ncols, "] ", x$data_type, " matrix\n", sep = "")
  cat("Number of observations: ", nrows, "\n", sep = "")
  cat("Number of time points: ", ncols, "\n", sep = "")
  type = switch(x$inter_type,
    "const_surv" = "Piece-wise Constant Survival",
    "linear_surv" = "Piece-wise Linear Survival",
    "const_haz" = "Piece-wise Constant Hazard",
  )
  cat("Interpolation type:", type)
}

#' @title Time Points
#' @description
#' Returns the time points associated with the object `x`.
#' @template param_x
#' @family methods
#' @export
times = function(x, ...) {
  UseMethod("times")
}

#' @rdname times
#' @return A `vector` of numeric values.
#' @export
times.survDistr = function(x, ...) {
  as.numeric(colnames(x$data))
}

#' @title Survival function
#' @description
#' Computes survival probabilities \eqn{S(t)} at the specified time points.
#' @details
#' If `x` is a [survDistr], and no time points are specified, returns the
#' stored survival `matrix`.
#' @template param_x
#' @family methods
#' @export
survival = function(x, ...) {
  UseMethod("survival")
}

#' @rdname survival
#' @template param_times
#' @return A `matrix` of survival probabilities.
#' @export
survival.survDistr = function(x, times = NULL) {
  if (x$data_type == "haz" || x$inter_type == "const_haz") {
    stop("Conversion from hazard to survival not yet implemented.")
  }

  if (is.null(times)) {
    if (x$data_type == "surv") {
      return(x$data) # survival matrix stored
    }
  }

  new_times = assert_numeric(times, lower = 0, any.missing = FALSE, min.len = 1)

  mat = rcpp_mat_interp(x = x$data, times = times(x), new_times = new_times, surv = TRUE,
                        constant = ifelse(x$inter_type == "const_surv", TRUE, FALSE))
  colnames(mat) = new_times

  mat
}

#' @title Cumulative Hazard function
#' @description
#' Computes the cumulative hazard at the specified time points as:
#' \eqn{H(t) = -log(S(t))}.
#' @details
#' If `x` is a [survDistr], and no time points are specified, the calculation of
#' \eqn{H(t)} uses the stored survival `matrix` \eqn{S(t)}.
#' @template param_x
#' @family methods
#' @export
cumhazard = function(x, ...) {
  UseMethod("cumhazard")
}

#' @rdname cumhazard
#' @template param_times
#' @template param_eps
#' @templateVar eps 1e-6
#' @return A `matrix` of cumulative hazards.
#' @export
cumhazard.survDistr = function(x, times = NULL, eps = 1e-6) {
  surv = survival(x, times = times)
  # Avoid S(t) = 0 by imputing a small epsilon to avoid log(0)
  surv[surv == 0] = eps

  -log(surv)
}

#' @title Hazard function
#' @description
#' Computes the hazard at the specified time points as: \eqn{h(t) = H(t) - H(t-1)}.
#' @details
#' If `x` is a [survDistr], and no time points are specified, the calculation of
#' \eqn{H(t)} uses the stored survival `matrix` \eqn{S(t)}.
#' @template param_x
#' @family methods
#' @export
hazard = function(x, ...) {
  UseMethod("hazard")
}

#' @rdname hazard
#' @template param_times
#' @template param_eps
#' @templateVar eps 1e-6
#' @return A `matrix` of hazards.
#' @export
hazard.survDistr = function(x, times = NULL, eps = 1e-6) {
  if (is.null(times)) {
    return(.rowwise_diffs(cumhazard(x, eps = eps)))
  }

  # Get the CDF matrix for the unique and sorted times
  utimes = sort(unique(times))
  haz = .rowwise_diffs(cumhazard(x, times = utimes)) # [nobs x utimes]

  # Create a mapping of `times` to `utimes`
  indx = match(times, utimes)

  # Extend hazard to match `times` using the mapping
  haz[, indx, drop = FALSE] # [nobs x times]
}

#' @title Cumulative Distribution Function
#' @description
#' Computes the cumulative distribution function \eqn{F(t)} at the specified time points.
#' \eqn{F(t)} is the probability that the event has occurred up until time \eqn{t}.
#' @details
#' If `x` is a [survDistr], and no time points are specified, returns
#' \eqn{F(t) = 1 - S(t)}, where \eqn{S(t)} is the stored survival matrix.
#' @template param_x
#' @family methods
#' @export
cdf = function(x, ...) {
  UseMethod("cdf")
}

#' @rdname cdf
#' @template param_times
#' @return A `matrix` of CDF values.
#' @export
cdf.survDistr = function(x, times = NULL) {
  if (is.null(times)) {
    return(1 - x$data)
  }

  new_times = assert_numeric(times, lower = 0, any.missing = FALSE, min.len = 1)
  mat = rcpp_mat_interp(x = 1 - x$data, times = times(x), new_times = new_times,
                        surv = FALSE, inter_type = x$inter_type)
  colnames(mat) = new_times

  mat
}

#' @title Probability Density Function
#' @description
#' Computes the probability density function \eqn{f(t)} at the specified time points.
#' \eqn{f(t)} is the probability of the event occurring at the specific time \eqn{t}.
#' @details
#' If `x` is a [survDistr], and no time points are specified, returns
#' \eqn{f(t) = F(t) - F(t-1)}, where \eqn{F(t)} is the cumulative distribution
#' function calculated directly from the stored survival matrix.
#' If new time points are specified, these are first ordered and duplicated values
#' are removed, so that the \eqn{F(t_{new})} can be calculated (always using linear
#' interpolation to avoid zero values for the pdf for time points lying in between
#' the ones we have S(t) in the stored survival matrix). pdf values for t_new < min(times)
#' or t_new > max(times) will be zero nonetheless.
#' @template param_x
#' @family methods
#' @export
pdf = function(x, ...) {
  UseMethod("pdf")
}

#' @rdname pdf
#' @template param_times
#' @return A `matrix` of PDF values.
#' @export
pdf.survDistr = function(x, times = NULL) {
  if (is.null(times)) {
    return(.rowwise_diffs(cdf(x)))
  }

  # Get the CDF matrix for the unique and sorted times
  utimes = sort(unique(times))
  pdf = .rowwise_diffs(cdf(x, times = utimes)) # [nobs x utimes]

  # Create a mapping of `times` to `utimes`
  indx = match(times, utimes)

  # Extend pdf to match `times` using the mapping
  pdf[, indx, drop = FALSE] # [nobs x times]
}
