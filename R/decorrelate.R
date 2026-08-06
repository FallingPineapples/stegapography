#' Flatten the y ~ x slope by duplicating influential points
#'
#' Drives the slope of `y ~ x` towards zero without moving any of the
#' existing points. On each iteration it finds the observation whose
#' influence on the slope (see [stats::influence()]) most closely cancels
#' the slope currently fit, and appends however many copies of that
#' observation are needed. It stops early once no positive number of
#' copies would help.
#'
#' Because rows are only ever appended, the result has at least as many
#' rows as the input, and usually a few more.
#'
#' @param data A data frame with numeric columns `x` and `y`.
#' @param max_iter Maximum number of duplication rounds.
#'
#' @return `data` with duplicated rows appended and `y` re-centred on zero.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(50), y = rnorm(50))
#' nrow(decorrelate_input(d))
decorrelate_input = function(data, max_iter=100) {
  data -> working_pineapple_data

  for (.i in 1:max_iter) {
    anti_model = lm(y ~ x, working_pineapple_data)
    gradient = -influence(anti_model)$coefficients[,2]
    point_id = which.min(abs(gradient - coef(anti_model)[[2]]))
    copies = round(coef(anti_model)[[2]] / gradient[[point_id]])
    if (copies <= 0) break
    working_pineapple_data %<>% .[c(1:nrow(.), rep(point_id, copies)),]
    # working_pineapple_data %<>% .[setdiff(1:nrow(.), point_id),]
  }

  working_pineapple_data %>% mutate(y = y - mean(y))
}

#' Impose a target correlation structure on a data frame
#'
#' Whitens every column, then re-imposes a correlation matrix that is the
#' identity except for the first two columns, whose observed 2x2 block is
#' preserved. Column means and standard deviations are restored afterwards.
#'
#' Method adapted from
#' <https://stackoverflow.com/questions/44930211/generate-uncorrelated-variables-each-well-correlated-with-existing-response-vari>.
#'
#' @param data A numeric data frame. The first two columns are the pair
#'   whose correlation is preserved; all other pairings are decorrelated.
#'
#' @return A tibble of the same dimensions as `data`.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(50), x0 = rnorm(50), x1 = rnorm(50))
#' round(cor(setcorrelate_output(d)), 3)
setcorrelate_output = function(data) {
  wpd = data
  mns = apply(wpd, 2, mean)
  sds = apply(wpd, 2, sd)

  wpd = apply(wpd, 2, scale)
  v.obs = cor(wpd)
  wpd = wpd %*% solve(chol(v.obs))
  want_cor = diag(nrow(v.obs))
  rownames(want_cor) = colnames(want_cor) = colnames(wpd)
  want_cor[1:2, 1:2] = v.obs[1:2, 1:2]
  wpd = wpd %*% chol(want_cor)

  wpd <- sweep(wpd, 2, sds, FUN="*")
  wpd <- sweep(wpd, 2, mns, FUN="+")

  as_tibble(wpd)
}
