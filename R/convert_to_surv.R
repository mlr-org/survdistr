#' @title Convert density/hazard to survival
#'
#' @description
#' Converts density or hazards from one of four input representations to survival
#' probabilities at the same anchor time points (no interpolation).
#'
#' @details
#' Let \eqn{t_1,\dots,t_B} denote the anchor time points,
#' \eqn{\Delta_j = t_j - t_{j-1}}, and \eqn{S_j = S(t_j)} the survival
#' probabilities at the anchors. The conversion depends on the value
#' of `input`:
#'
#' * `"disc_dens"`: \eqn{S_j = 1 - \sum_{k=1}^j \tilde f_k}
#'
#' * `"disc_haz"`: \eqn{S_j = \prod_{k=1}^j (1 - \tilde h_k)}
#'
#' * `"cont_dens"`: \eqn{S_j = 1 - \sum_{k=1}^j f_k \Delta_k}
#'
#' * `"cont_haz"`: \eqn{S_j = \exp\!\left(-\sum_{k=1}^j \lambda_k \Delta_k\right)}
#'
#' @section Validation:
#' If `check = TRUE`, we validate that the input is a proper discrete density/hazard matrix
#' or vector using [assert_prob()].
#' For continuous hazards/densities, we only check that the input is a non-negative numeric
#' matrix/vector.
#'
#' @param x (`numeric()` | `matrix()`)\cr
#'  Input vector or matrix (rows = observations, columns = time points).
#' @param times (`numeric()` | `NULL`)\cr
#'  Anchor time points. If `NULL`, extracted from names/colnames of `x`.
#' @param input (`character(1)`)\cr
#'  Input type. One of `"disc_haz"`, `"disc_dens"`, `"cont_haz"` or `"cont_dens"`.
#' @param check (`logical(1)`)\cr
#'  If `TRUE` (default), run \emph{input} validation checks.
#'  Disable only if you know the input is valid and want to skip checks for speed.
#' @param clamp_surv (`logical(1)`)\cr
#'  If `TRUE`, clamp survival probabilities to `[eps, 1]` to avoid numerical issues.
#' @param eps (`numeric(1)`)\cr
#'  Small value used to clamp near-zero survival probabilities if `clamp_surv = TRUE`.
#'
#' @return A numeric vector or matrix of survival probabilities with the same
#'  dimensions as `x`.
#' @examples
#' # Continuous hazard => survival
#' haz_cont = c(0.02, 0.1, 0.2, 0.15)
#' times = c(0, 1, 2, 3)
#' convert_to_surv(haz_cont, times = times, input = "cont_haz")
#'
#' # Discrete hazard => survival
#' haz_disc = c(0.1, 0.2, 0.15)
#' times = c(1, 2, 3)
#' convert_to_surv(haz_disc, times = times, input = "disc_haz")
#'
#' @export
convert_to_surv = function(x, times = NULL, input = "cont_haz", check = TRUE,
                           clamp_surv = FALSE, eps = 1e-6) {
  check = assert_flag(check)
  clamp_surv = assert_flag(clamp_surv)
  input = assert_choice(input, c("disc_dens", "disc_haz", "cont_dens", "cont_haz"))
  times = extract_times(x, times)
  is_mat = is.matrix(x)
  x_mat = if (is_mat) x else matrix(x, nrow = 1)

  if (check) {
    if (startsWith(input, "disc")) {
      assert_prob(x = x_mat, times = times, type = if (input == "disc_dens") "dens" else "haz")
    } else {
      assert_matrix(x_mat, mode = "numeric", any.missing = FALSE, min.rows = 1, min.cols = 1)
      if (any(x_mat < 0)) {
        stop("For continuous hazards/densities, all values must be non-negative.")
      }
    }
  }

  surv = switch(input,
    "disc_dens" = c_disc_dens_to_surv_mat(x_mat),
    "disc_haz"  = c_disc_haz_to_surv_mat(x_mat),
    "cont_dens" = c_cont_dens_to_surv_mat(x_mat, times),
    "cont_haz"  = c_cont_haz_to_surv_mat(x_mat, times)
  )

  if (clamp_surv) {
    surv = c_clamp_surv(surv, eps = eps)
  }

  if (!is_mat) {
    surv = surv[1, ]
  }

  surv
}
