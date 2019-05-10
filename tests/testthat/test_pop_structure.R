context("Population Structure")
library(dplyr)

test_that("the population structure in pop and other variables match", {
  expect_true(all(pop[pop1,'i'] == 'black')) # test subpop dimension
  expect_true(all(pop[pop2,'i'] == 'white'))
  expect_true(all(pop[pop3,'i'] == 'hispanic'))
  expect_true(all(pop[pop4,'i'] == 'msm-hivneg'))
  expect_true(all(pop[pop5,'i'] == 'msm-hivpos'))

  expect_true(all(pop[males,'j'] == 'male')) # test sex dimension
  expect_true(all(pop[females,'j'] == 'female'))

  expect_true(all(pop[y.m, 'l'] == 'young')) # test age dimension
  expect_true(all(pop[y.f, 'l'] == 'young')) # test age dimension
  expect_true(all(pop[o.m, 'l'] == 'old')) # test age dimension
  expect_true(all(pop[o.m, 'l'] == 'old')) # test age dimension

	expect_true(all(pop[pop$k == 'low', 'index'] %% 2 == 1)) # low activity are odd
	expect_true(all(pop[pop$k == 'high', 'index'] %% 2 == 0)) # high activity are even

})
