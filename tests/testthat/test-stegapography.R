picture = function(n = 50) {
  data.frame(x = rnorm(n), y = rnorm(n))
}

test_that("output has n_x + 1 columns, named y and x0..x{n_x-1}", {
  set.seed(1)
  out = stegapography(picture(), n_x = 5)

  expect_named(out, c("y", "x0", "x1", "x2", "x3", "x4"))
})

test_that("n_x is respected", {
  set.seed(2)
  expect_length(stegapography(picture(), n_x = 3), 4)
  expect_length(stegapography(picture(), n_x = 8), 9)
})

test_that("rows are appended, never dropped", {
  set.seed(3)
  d = picture(50)

  expect_gte(nrow(stegapography(d)), nrow(d))
  expect_gte(nrow(decorrelate_input(d)), nrow(d))
})

test_that("decorrelate_input flattens the y ~ x slope", {
  set.seed(4)
  d = picture(50)

  before = abs(coef(lm(y ~ x, d))[[2]])
  after = abs(coef(lm(y ~ x, decorrelate_input(d)))[[2]])

  expect_lt(after, before)
})

test_that("setcorrelate_output decorrelates everything but the first pair", {
  set.seed(5)
  d = data.frame(y = rnorm(50), x0 = rnorm(50), x1 = rnorm(50), x2 = rnorm(50))
  v = cor(setcorrelate_output(d))

  # the y/x0 block is preserved
  expect_equal(v[1:2, 1:2], cor(d)[1:2, 1:2])
  # every other off-diagonal correlation is gone
  expect_equal(v[3:4, 3:4], diag(2), ignore_attr = TRUE)
  expect_equal(v[1:2, 3:4], matrix(0, 2, 2), ignore_attr = TRUE)
})

test_that("mangle keeps only the columns it created", {
  out = mangle(data.frame(x = 1:3, y = 4:6), b = y, a = x)

  expect_named(out, c("b", "a"))
  expect_equal(out$b, 4:6)
})

test_that("v2v_rotation_matrix rotates x onto y", {
  from = c(1, 0, 0)
  to = normalize(c(1, 2, 3))

  expect_equal(drop(from %*% v2v_rotation_matrix(from, to)), to)
  expect_error(v2v_rotation_matrix(c(1, 1, 1), to), "is_normalized")
})

test_that("normalize gives unit length", {
  expect_equal(normalize(c(3, 4)), c(0.6, 0.8))
  expect_true(is_normalized(normalize(rnorm(10))))
  expect_false(is_normalized(c(3, 4)))
})
