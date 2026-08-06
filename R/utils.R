#' Mutate and keep only the new columns
#'
#' A [dplyr::mutate()] that discards every column it was not asked to
#' create, so `mangle(d, y = y, x0 = x)` returns a two-column frame named
#' `y` and `x0`.
#'
#' This uses non-standard evaluation to recover the names of `...`, so it
#' cannot be called programmatically with a pre-built list of expressions.
#'
#' @param .data A data frame.
#' @param ... Named expressions, as passed to [dplyr::mutate()]. Only the
#'   columns named here survive.
#'
#' @return A data frame containing only the columns named in `...`, in that
#'   order.
#' @export
#' @examples
#' mangle(data.frame(x = 1:3, y = 4:6), doubled = x * 2)
mangle = function(.data, ...) {
  nms = names(substitute(stop(...)))[-1]

  mutate(.data, ...) %>%
    select(one_of(nms))
}

#' Scale a vector to unit length
#'
#' @param x A numeric vector.
#'
#' @return `normalize()` returns `x` divided by its Euclidean norm.
#'   `is_normalized()` returns a single `TRUE` or `FALSE`, testing unit
#'   length up to [all.equal()] tolerance.
#' @export
#' @examples
#' normalize(c(3, 4))
#' is_normalized(normalize(c(3, 4)))
normalize = function(x) {
  x / sqrt(drop(x %*% x))
}

#' @rdname normalize
#' @export
is_normalized = function(x) {
  isTRUE(all.equal(drop(x %*% x), 1))
}
