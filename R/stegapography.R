mangle = function(.data, ...) {
  nms = names(substitute(stop(...)))[-1]

  mutate(.data, ...) %>%
    select(one_of(nms))
}

stegapograpy = function(dithered_pineapple_data, n_x=4) {
  dithered_pineapple_data %>%
  	mutate(across(everything(), scale)) -> rescaled_pineapple_data

  rescaled_pineapple_data %>% decorrelate_input() -> display_pineapple_data

  display_pineapple_data %>% mutate(y = x + y) -> data_pineapple_data

  setNames(nm=paste0("x", 1:(n_x-1))) %>%
    map(\() rnorm(nrow(data_pineapple_data))) -> random_data

  data_pineapple_data %>%
    mangle(y = y, x0 = x) %>%
    bind_cols(random_data) ->
    expanded_pineapple_data

  expanded_pineapple_data %>% setcorrelate_output() -> decor_pineapple_data

  x0_basis = c(0, 1, rep.int(0, n_x-1))
  target = c(0, sample(c(1, -1), n_x, replace=TRUE))
  final_transform = v2v_rotation_matrix(x0_basis, target)

  (as.matrix(decor_pineapple_data) %*% final_transform) %>% as.data.frame()
}

decorrelate_input = function(rescaled_pineapple_data) {
  rescaled_pineapple_data -> working_pineapple_data

  for (i in 1:100) {
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

# reduce more of the correlations
# https://stackoverflow.com/questions/44930211/generate-uncorrelated-variables-each-well-correlated-with-existing-response-vari
setcorrelate_output = function(expanded_pineapple_data) {
  wpd = expanded_pineapple_data
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

# https://en.wikipedia.org/wiki/Rotation_matrix#Vector_to_vector_formulation
v2v_rotation_matrix = function(x, y) {
  I = diag(nrow=length(x))
  m = (x %o% y) - (y %o% x)
  I + m + (1/(1+(x %*% y)))*(m %*% m)
}
