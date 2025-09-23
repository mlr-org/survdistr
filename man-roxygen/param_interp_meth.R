#' @param interp_meth (`character(1)`)\cr
#'  Interpolation method to use when requesting the quantity of interest at time
#'  points different than the ones in the stored object (accessible via the `times` method).
#'  Currently supported interpolation methods include `"const_surv"` (default),
#'  `"linear_surv"` and `"const_haz"`. See details.
