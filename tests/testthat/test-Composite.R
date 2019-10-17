context("Generic tests")
test_that("Sanity checks", {

  expect_error(
    Composite(var = array(1:20, dim = c(2, 5, 2)), c(1, 1, 0)),
    "Temporal dimension of var is not equal to length of occ.")

  expect_warning(
    Composite(var = array(1:40, dim = c(2, 5, 4)), c(1, 2, 2, 2)),
    "Composite 1 has length 1 and pvalue is NA.")

  var <- array(rep(c(1, 3, 2, 1, 2), 8), dim = c(x = 2, y = 4, time = 5))
  occ <- c(1, 2, 2, 2, 1)
  output <- c(x = 2, y = 4, 2)  #dim(asd$composite)
  expect_equal(
    dim(Composite(var, occ)$composite),
    output
  )
  output <- c(1.5, 2.0, 2.5, 2.0)
  expect_equal(
    Composite(var, occ)$composite[1, , 1],
    output
  )

  var <- array(rep(c(1, 3, 2, 1, 2, 2), 8), dim = c(x = 2, y = 4, time = 6))
  occ <- c(1, 1, 2, 2, 3, 3)
  output <- matrix(c(1.5, 2.5, 1.5, 2.0, 2.0, 1.5, 1.5, 2.5), 2, 4)
  expect_equivalent(
    Composite(var, occ)$composite[, , 2],
    output
  )

})

