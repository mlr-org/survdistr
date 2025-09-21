#' @title Survival Distribution Container
#' @name survDistr
#'
#' @description
#' [survDistr] is an [R6][R6::R6Class] specialized container designed for storing
#' and managing prediction outputs from survival models in single-event settings.
#' This includes models such as Cox proportional hazards, random survival forests,
#' and other classical or machine learning-based survival estimators.
#'
#' The main prediction data type can be a survival or a hazard matrix, where
#' **rows represent observations and columns represent time points**.
#'
#' @template param_times
#' @template param_add_times
#' @template param_eps
#' @templateVar eps 1e-6
#'
#' @details
#' The input matrix (survival probabilities \eqn{S(t)} or hazard \eqn{h(t)})
#' is stored internally and accessed by the `$data` field.
#' the interpolation type needed for the
#' public methods is stored in the `$interp_meth` slot.
#'
#' During construction, the function [assert_prop_matrix()] is used to validate
#' the input data matrix according to the given `data_type`.
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
#' x$interp_meth
#'
#' # time points
#' x$times
#'
#' # S(t) at given time points (constant interpolation)
#' x$survival(times = c(10, 30, 42, 50))
#' # same but with linear interpolation
#' x$interp_meth = "linear_surv"
#' x$survival(times = c(10, 30, 42, 50))
#' # time points can be unordered and duplicated
#' x$survival(times = c(10, 30, 10, 50))
#'
#' # Cumulative hazard
#' x$cumhazard()
#'
#' @export
survDistr = R6Class(
  "survDistr",

  public = list(
    #' @field times (`numeric`])\cr
    #'  Numeric vector of time points corresponding to columns of `data`.
    times = NULL,
    #' @field data_type (`character(1)`)\cr
    #'  Either `"surv"` for survival or `"haz"` for hazard matrices.
    data_type = NULL,
    #' @field interp_meth (`character(1)`)\cr
    #'  Interpolation method; one of `"const_surv"`, `"linear_surv"`, or `"const_haz"`.
    interp_meth = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param x (`matrix`)\cr
    #'  A numeric matrix of either survival probabilities (values between 0 and 1) or
    #'  hazard values (non-negative values).
    #'  Column names must correspond to time points if `times` is `NULL`.
    #' @template param_data_type
    #' @template param_interp_meth
    #' @param times (`numeric(1)`)\cr Numeric vector of time points for matrix `x`,
    #'  must match the number of columns.
    #' @param ... currently not used
    initialize = function(x, times = NULL, data_type = "surv", interp_meth = "const_surv", ...) {
      # Validate input type
      assert_choice(data_type, c("surv", "haz"))

      if (data_type == "surv") {
        times = assert_prob_matrix(x, times, type = "surv")
      } else {
        stop("Input hazard matrix not yet supported.")
      }

      dimnames(x) = NULL # no need to keep these
      private$.mat = x # store data matrix

      # Validate interpolation method
      assert_choice(interp_meth, c("const_surv", "linear_surv", "const_haz"))
      if (interp_meth == "const_haz") {
        stop("Constant hazard interpolation not yet implemented.")
      }

      # TODO: check is this is needed at all
      if (data_type == "haz" && interp_meth == "linear_surv") {
        stop("Hazard data and piece-wise linear survival interpolation is not supported.")
      }

      # Fill in public fields
      self$times = times
      self$data_type = data_type
      self$interp_meth = interp_meth
    },

    #' @description
    #' Displays summary information about a [survDistr] object, including
    #' the number of observations and time points.
    print = function() {
      nrows = nrow(private$.mat)
      ncols = ncol(private$.mat)

      data_type = switch(self$data_type,
                         "surv" = "survival",
                         "haz" = "hazard")
      cat("A [", nrows, " x ", ncols, "] ", data_type, " matrix\n", sep = "")
      cat("Number of observations: ", nrows, "\n", sep = "")
      cat("Number of time points: ", ncols, "\n", sep = "")
      interp_meth = switch(self$interp_meth,
                           "const_surv" = "Piece-wise Constant Survival",
                           "linear_surv" = "Piece-wise Linear Survival",
                           "const_haz" = "Piece-wise Constant Hazard")
      cat("Interpolation method:", interp_meth, "\n")
      invisible(self)
    },

    #' @description
    #' Return the stored data matrix.
    #' @return (`matrix`)
    data = function(add_times = TRUE) {
      assert_flag(add_times)

      mat = private$.mat
      if (add_times) {
        colnames(mat) = as.character(self$times)
      }

      mat
    },

    #' @description
    #' Computes survival probabilities \eqn{S(t)} at the specified time points.
    #' Uses [mat_interp()].
    #'
    #' @return a `matrix` of survival probabilities
    survival = function(times = NULL, add_times = TRUE) {
      mat_interp(
        x = private$.mat,
        times = self$times,
        eval_times = times,
        constant = self$interp_meth == "const_surv",
        type = "surv",
        add_times = add_times,
        check = FALSE # input `x` is already checked in initialize()
      )
    },

    #' @description
    #' Computes the cumulative distribution function \eqn{F(t) = 1 - S(t)} at the specified time points.
    #' \eqn{F(t)} is the probability that the event has occurred up until time \eqn{t}.
    #' Uses [mat_interp()].
    #'
    #' @return a cdf `matrix`.
    cdf = function(times = NULL, add_times = TRUE) {
      mat_interp(
        x = 1 - private$.mat, # convert survival => CDF
        times = self$times,
        eval_times = times,
        constant = self$interp_meth == "const_surv",
        type = "cdf",
        add_times = add_times,
        check = FALSE # input `x` is already checked in initialize()
      )
    },

    #' @description
    #' Computes the cumulative hazard at the specified time points as:
    #' \eqn{H(t) = -log(S(t))}.
    #'
    #' @return a `matrix` of cumulative hazards.
    cumhazard = function(times = NULL, add_times = TRUE, eps = 1e-6) {
      surv_mat = self$survival(times = times, add_times = add_times)
      surv_mat[surv_mat == 0] = eps
      -log(surv_mat)
    },

    #' @description
    #' Computes the hazard at the specified time points as: \eqn{h(t) = H(t) - H(t-1)}.
    #'
    #' @return a hazard `matrix`.
    hazard = function(times = NULL, eps = 1e-6) {
      if (is.null(times)) {
        return(rowwise_diffs(self$cumhazard(eps = eps)))
      }

      utimes = sort(unique(times))
      haz = rowwise_diffs(self$cumhazard(times = utimes))

      indx = match(times, utimes)
      haz[, indx, drop = FALSE]
    },

    #' @description
    #' Computes the probability density function \eqn{f(t)} at the specified time points.
    #' \eqn{f(t)} is the probability of the event occurring at the specific time \eqn{t}.
    #' For constant survival interpolation, \eqn{f(t) = F(t) - F(t-1)}, where
    #' \eqn{F(t)} is the cumulative distribution.
    #'
    #' @return a pdf `matrix`.
    pdf = function(times = NULL) {
      if (self$interp_meth != "const_surv") {
        stop("Only implemented for constant survival interpolation.")
      }

      if (is.null(times)) {
        return(rowwise_diffs(self$cdf()))
      }

      utimes = sort(unique(times))
      pdf_mat = rowwise_diffs(self$cdf(times = utimes))

      indx = match(times, utimes)
      pdf_mat[, indx, drop = FALSE]
    }
  ),

  private = list(
    .mat = NULL
  )
)
