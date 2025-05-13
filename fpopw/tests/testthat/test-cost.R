library(testthat)
test_that("default is no model", {
  data.vec <- 1:5
  penalty <- 10
  fit <- fpopw::Fpop(data.vec, penalty)
  cum.squares <- cumsum(data.vec^2)
  expected.cost <- sapply(
    1:length(data.vec),
    function(i)sum((data.vec[1:i]-mean(data.vec[1:i]))^2))
  computed.cost <- fit$cost+cum.squares-penalty
  expect_equal(computed.cost, expected.cost)
  expect_identical(fit$model, NULL)
  ##A$J.est <- A$cost[n] - A$K*lambda + sum(x^2)
  expect_equal(fit$J.est, sum((mean(data.vec)-data.vec)^2))
})
test_that("verbose_file yields model", {
  fit <- fpopw::Fpop(1:5, 10, verbose_file=tempfile())
  expect_is(fit$model, "data.table")
  expect_is(fit$model$cost, "numeric")
  expect_is(fit$model$change_i, "integer")
})
