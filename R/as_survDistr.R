#' @title Coerce Object to [survDistr]
#'
#' @description
#' S3 generic to coerce supported objects into a [survDistr] object.
#'
#' @param x Object to convert.
#' @param ... Additional arguments passed to methods.
#' @return A [survDistr] object.
#' @export
as_survDistr = function(x, ...) {
  UseMethod("as_survDistr")
}

#' @rdname as_survDistr
#' @title Convert a matrix of survival probabilities into a [survDistr] object.
#' @param times (`numeric`)
#'   Numeric vector of time points corresponding to columns of `x`.
#'   If `NULL`, column names of `x` are used.
#' @param method (`character(1)`)\cr
#'   Interpolation method passed to [survDistr] constructor.
#' @param check (`logical(1)`)\cr
#'   Whether to validate `x` and `times`.
#' @param trim_duplicates (`logical(1)`)\cr
#'  Whether to remove duplicate S(t) values and corresponding time points.
#' @export
as_survDistr.matrix = function(x, times = NULL, method = "const_surv", check = TRUE,
                               trim_duplicates = FALSE) {
  survDistr$new(x = x, times = times, method = method, check = check, trim_duplicates = trim_duplicates)
}

#' @rdname as_survDistr
#' @title Return [survDistr] objects unchanged.
#' @export
as_survDistr.survDistr = function(x, ...) {
  x
}

#' @rdname as_survDistr
#' @export
as_survDistr.default = function(x, ...) {
  stop(
    "No as_survDistr() method for objects of class ",
    paste(class(x), collapse = "/"),
    ".",
    call. = FALSE
  )
}
