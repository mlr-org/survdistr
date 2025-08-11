#' @param inter_type (`character(1)`)\cr
#' The type of interpolation to use for time points different than the ones in
#' the stored object (accessible via the [times()] method).
#' Currently supported interpolation methods include `"const_surv"` (default),
#' `"linear_surv"` and `"const_haz"`. See details.
#'
