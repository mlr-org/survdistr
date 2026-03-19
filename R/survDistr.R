#' @title Survival Distribution Container
#' @name survDistr
#'
#' @description
#' [survDistr] is an [R6][R6::R6Class] specialized container designed for storing
#' and managing prediction outputs from survival models in single-event settings.
#' This includes models such as Cox proportional hazards, random survival forests,
#' and other classical or machine learning-based survival estimators.
#'
#' The main prediction data type is survival matrix, where
#' **rows represent observations and columns represent time points**.
#'
#' @template param_times
#' @template param_add_times
#' @template param_rows
#' @template param_eps
#'

#' @details
#' Key design features:
#' - The survival matrix is stored internally and can be accessed using the `$data()` method.
#' - The `$times` active field provides the time points corresponding to the matrix columns.
#' - The interpolation method is controlled via the `$method` active field.
#' - Survival-related quantities (e.g., distribution, density, hazard functions) are
#' interpolated using the [interp()] function.
#' - The [assert_prob()] function validates the input data matrix during construction if
#' `check` is set to `TRUE`.
#' - Use the `$filter()` method to subset observations in-place by filtering rows of the
#' stored matrix.
#' - Use `trim_duplicates = TRUE` in the constructor to remove flat survival segments (repeated
#' values across time points) with a pre-specified tolerance (for a more controlled
#' pre-processing, see [trim_duplicates()]).
#'
#' @examples
#' # generate survival matrix
#' mat = matrix(data = c(1,0.6,0.4,0.8,0.8,0.7), nrow = 2,
#'              ncol = 3, byrow = TRUE)
#' times = c(12, 34, 42)
#' x = survDistr$new(mat, times)
#' x
#'
#' # stored survival matrix
#' x$data()
#'
#' # interpolation method
#' x$method
#'
#' # time points
#' x$times
#'
#' eval_times = c(10, 30, 42, 50)
#' # S(t) at given time points (constant interpolation)
#' x$survival(times = eval_times)
#' # same but with linear interpolation
#' x$method = "linear_surv"
#' x$survival(times = eval_times)
#'
#' # Cumulative hazard
#' x$cumhazard(times = eval_times)
#'
#' # Density
#' x$density(times = eval_times)
#'
#' # Hazard
#' x$hazard(times = eval_times)
#'
#' @export
survDistr = R6Class(
  "survDistr",

  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param x (`matrix`)\cr
    #'  A numeric matrix of survival probabilities (values between 0 and 1).
    #'  Column names must correspond to time points if `times` is `NULL`.
    #' @param times (`numeric(1)`)\cr Numeric vector of time points for matrix `x`,
    #'  must match the number of columns.
    #' @template param_method
    #' @template param_check
    #' @template param_trim_duplicates
    initialize = function(x, times = NULL, method = "const_surv", check = TRUE,
                          trim_duplicates = FALSE) {
      assert_flag(check)
      assert_flag(trim_duplicates)
      method = map_interp_method(method) # const_* aliases
      private$.method = method

      # remove flat S(t) segments
      if (trim_duplicates) {
        trimmed = trim_duplicates(x, times = times)
        x = trimmed$x
        times = trimmed$times
      }

      if (check) {
        times = assert_prob_matrix(x, times, type = "surv")
      } else {
        times = extract_times(x, times)
      }
      private$.times = times

      dimnames(x) = NULL # no need to keep these
      private$.mat = x # store data matrix
    },

    #' @description
    #' Displays summary information about a [survDistr] object, including
    #' the number of observations and time points.
    print = function() {
      nrows = nrow(private$.mat)
      ncols = ncol(private$.mat)

      cat("A [", nrows, " x ", ncols, "] survival matrix\n", sep = "")
      cat("Number of observations: ", nrows, "\n", sep = "")
      cat("Number of time points: ", ncols, "\n", sep = "")
      method = switch(
        self$method,
        "const_surv" = "Piecewise Constant Survival",
        "const_dens" = "Piecewise Constant Density (Linear Survival)",
        "const_haz"  = "Piecewise Constant Hazard (Exponential Survival)"
      )
      cat("Interpolation method:", method, "\n")
      invisible(self)
    },

    #' @description
    #' Return the stored data matrix.
    #'
    #' @return (`matrix`)
    data = function(rows = NULL, add_times = TRUE) {
      assert_flag(add_times)

      mat = private$.filter_mat(rows)
      if (add_times) {
        colnames(mat) = as.character(self$times)
      }

      mat
    },

    #' @description
    #' Filters observations **in-place** by subsetting rows of the stored matrix.
    #'
    #' @return Invisibly returns the `survDistr` object itself.
    filter = function(rows = NULL) {
      if (is.null(rows)) {
        return(invisible(self))
      }

      private$.mat = private$.filter_mat(rows)
      invisible(self)
    },

    #' @description
    #' Computes survival probabilities \eqn{S(t)} at the specified time points.
    #'
    #' @return a `matrix` of survival probabilities (rows = observations, columns = time points).
    survival = function(rows = NULL, times = NULL, add_times = TRUE) {
      interp(
        x = private$.filter_mat(rows),
        times = self$times,
        eval_times = times,
        method = self$method,
        output = "surv",
        add_times = add_times,
        check = FALSE # input `x` is already checked in initialize()
      )
    },

    #' @description
    #' Computes the cumulative distribution function \eqn{F(t) = 1 - S(t)} or CDF at the specified time points.
    #' \eqn{F(t)} is the probability that the event has occurred up until time \eqn{t}.
    #'
    #' @return a `matrix` of CDF values (rows = observations, columns = time points).
    distribution = function(rows = NULL, times = NULL, add_times = TRUE) {
      interp(
        x = private$.filter_mat(rows),
        times = self$times,
        eval_times = times,
        method = self$method,
        output = "cdf",
        add_times = add_times,
        check = FALSE # input `x` is already checked in initialize()
      )
    },

    #' @description
    #' Computes the probability density function \eqn{f(t)} or PDF at the specified time points.
    #' \eqn{f(t) = \frac{d}{dt} F(t)} is the probability of the event occurring at the specific
    #' time \eqn{t}.
    #'
    #' @return a `matrix` of PDF values (rows = observations, columns = time points).
    density = function(rows = NULL, times = NULL, add_times = TRUE) {
      interp(
        x = private$.filter_mat(rows),
        times = self$times,
        eval_times = times,
        method = self$method,
        output = "density",
        add_times = add_times,
        check = FALSE # input `x` is already checked in initialize()
      )
    },

    #' @description
    #' Computes the cumulative hazard (accumulated risk up to time \eqn{t}) at the specified time
    #' points as \eqn{H(t) = -log(S(t))}.
    #'
    #' @return a `matrix` of cumulative hazards (rows = observations, columns = time points).
    cumhazard = function(rows = NULL, times = NULL, add_times = TRUE, eps = 1e-12) {
     interp(
        x = private$.filter_mat(rows),
        times = self$times,
        eval_times = times,
        method = self$method,
        output = "cumhaz",
        add_times = add_times,
        check = FALSE, # input `x` is already checked in initialize()
        eps = eps
      )
    },

    #' @description
    #' Computes the hazard \eqn{h(t) = \frac{f(t)}{S(t)}} at the specified time points.
    #' Hazard is the conditional instantaneous event rate at time \eqn{t} given
    #' survival up to time \eqn{t}.
    #'
    #' @return a `matrix` of hazard values (rows = observations, columns = time points).
    hazard = function(rows = NULL, times = NULL, add_times = TRUE) {
      interp(
        x = private$.filter_mat(rows),
        times = self$times,
        eval_times = times,
        method = self$method,
        output = "hazard",
        add_times = add_times,
        check = FALSE # input `x` is already checked in initialize()
      )
    }
  ),

  active = list(
    #' @field times (`numeric`)\cr
    #'  Numeric vector of time points corresponding to columns of `data`. Read-only.
    times = function(rhs) {
      if (missing(rhs)) return(private$.times)
      stop("`times` is read-only.")
    },

    #' @field method (`character(1)`)\cr
    #'  Interpolation method; one of `"const_surv"` (default), `"const_dens"` (alias: `"linear_surv"`)
    #'  and `"const_haz"` (alias: `"exp_surv"`).
    method = function(rhs) {
      if (missing(rhs)) return(private$.method)
      private$.method = map_interp_method(rhs)
    }
  ),

  private = list(
    .mat = NULL,
    .times = NULL,
    .method = NULL,
    .filter_mat = function(rows = NULL) {
      # check rows and return filtered matrix
      if (is.null(rows)) {
        return(private$.mat)
      }

      n = nrow(private$.mat)
      if (is.logical(rows)) {
        rows = assert_logical(rows, any.missing = FALSE, len = n)
      } else {
        rows = assert_integerish(rows, lower = 1L, upper = n, unique = TRUE, sorted = TRUE, any.missing = FALSE)
      }

      private$.mat[rows, , drop = FALSE]
    }
  )
)
