#' Hide a picture in a regression's residuals
#'
#' Takes a two-column picture and expands it into `n_x + 1` columns that
#' look like unremarkable predictors: the point cloud is decorrelated, the
#' picture is folded into `y`, noise columns are added, and the whole thing
#' is rotated so the signal is smeared across every `x` column rather than
#' sitting in one of them.
#'
#' The picture is recovered by regressing `y` on the `x` columns and
#' plotting the residuals.
#'
#' @param data A data frame with numeric columns `x` and `y`, holding the
#'   picture to hide.
#' @param n_x Number of `x` columns to produce.
#' @param max_iter Passed to [decorrelate_input()].
#'
#' @return A data frame with `n_x + 1` columns named `y`, `x0`, `x1`, ...,
#'   `x{n_x - 1}`. It has *more rows than `data`*, because
#'   [decorrelate_input()] flattens the input slope by duplicating points.
#' @export
#' @examples
#' set.seed(1)
#' picture <- data.frame(x = rnorm(50), y = rnorm(50))
#' hidden <- stegapography(picture)
#' round(cor(hidden), 3)
#' plot(resid(lm(y ~ ., hidden)) ~ fitted(lm(y ~ ., hidden)))
stegapography = function(data, n_x=5, max_iter=100) {
  data %>% mutate(across(everything(), scale)) -> rescaled_pineapple_data

  rescaled_pineapple_data %>% decorrelate_input(max_iter) -> display_pineapple_data

  display_pineapple_data %>% mutate(y = x + y) -> data_pineapple_data

  setNames(nm=paste0("x", 1:(n_x-1))) %>%
    map(\(.x) rnorm(nrow(data_pineapple_data))) -> random_data

  data_pineapple_data %>%
    mangle(y = y, x0 = x) %>%
    bind_cols(random_data) ->
    expanded_pineapple_data

  expanded_pineapple_data %>% setcorrelate_output() -> decor_pineapple_data

  x0_basis = c(0, 1, rep.int(0, n_x-1))
  # TODO: target is obvious in correlations, remove oddities
  target = c(0, sample(c(1, -1), n_x, replace=TRUE)) %>% normalize()
  final_transform = v2v_rotation_matrix(x0_basis, target)
  dimnames(final_transform) <- list(names(decor_pineapple_data), names(decor_pineapple_data))

  (as.matrix(decor_pineapple_data) %*% final_transform) %>% as.data.frame()
}
