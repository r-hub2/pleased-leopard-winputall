test_that("compare rp_model_fun", {
  expect_equal(rp_model_fun( beta=c(1,2,3), opt=list(distrib="lognormal")), as.matrix(exp(c(1,2,3))))
})
