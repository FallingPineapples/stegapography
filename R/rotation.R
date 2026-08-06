#' Rotation matrix taking one unit vector onto another
#'
#' Vector-to-vector formulation, see
#' <https://en.wikipedia.org/wiki/Rotation_matrix#Vector_to_vector_formulation>.
#'
#' @param x,y Unit vectors of the same length. Both are checked with
#'   [is_normalized()].
#'
#' @return A square rotation matrix `R` of side `length(x)` such that
#'   `drop(x %*% R)` equals `y`.
#' @export
#' @examples
#' R <- v2v_rotation_matrix(c(1, 0, 0), normalize(c(1, 1, 1)))
#' drop(c(1, 0, 0) %*% R)
v2v_rotation_matrix = function(x, y) {
  stopifnot(length(x) == length(y))
  stopifnot(is_normalized(x))
  stopifnot(is_normalized(y))
  I = diag(nrow=length(x))
  m = (x %o% y) - (y %o% x)
  I + m + (1/(1+drop(x %*% y)))*(m %*% m)
}
