#' stegapography: Hide Images in Residual Plots
#'
#' Takes a two-column picture and disguises it as an unremarkable-looking
#' set of predictors. The picture is invisible in the raw data and in the
#' pairwise correlations, but reappears in the residuals of a linear model
#' fit to the result.
#'
#' The entry point is [stegapography()]; the pipeline stages it calls are
#' documented individually.
#'
#' @importFrom dplyr mutate select one_of across everything bind_cols as_tibble
#' @importFrom purrr map
#' @importFrom magrittr %>% %<>%
#' @importFrom stats lm influence coef rnorm sd cor setNames
#' @keywords internal
"_PACKAGE"
