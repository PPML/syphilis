context("The package should automatically provide optimized parameter vectors.")

state <<- 'LA'
load.start()
d <- dLogPosterior(theta_la)
criterion <- function(d) is.numeric(d) && is.finite(d) && d > -1e32

test_that("dLogPosterior(theta_la) yields a finite number > -1e32.", {
  expect_true(criterion(d))
})

state <<- 'MA'
load.start()
d <- dLogPosterior(theta_ma)

test_that("dLogPosterior(theta_ma) yields a finite number > -1e32.", {
  expect_true(criterion(d))
})
