test_that("likelihood ratios use log-likelihoods consistently for either sign", {
  for (pair in list(c(-6, -8), c(8, 6))) {
    for (df in c(1, 3)) {
      result <- pbiasr:::.ak_likelihood_ratio(pair[1], pair[2], df)
      expect_equal(result$statistic, 4)
      expect_equal(result$p_value, pchisq(4, df, lower.tail = FALSE))
    }
  }
  expect_equal(pbiasr:::.ak_likelihood_ratio(10, 10, 1)$p_value, 1)
  expect_equal(pbiasr:::.ak_likelihood_ratio(10 - 1e-10, 10, 1)$statistic, 0)
  expect_equal(pbiasr:::.ak_likelihood_ratio(8, 10, 1)$status, "failed")
  expect_true(is.na(pbiasr:::.ak_likelihood_ratio(Inf, 10, 1)$p_value))
})

test_that("unit publication weights recover the ordinary normal likelihood", {
  d <- data.frame(id = 1:6, effect = c(-0.7, -0.2, -0.05, 0.1, 0.3, 0.8),
                  se = c(0.1, 0.2, 0.15, 0.1, 0.2, 0.15), constant = 1)
  expected <- dnorm(d$effect, mean = 0.05, sd = sqrt(d$se^2 + 0.04), log = TRUE)
  expect_equal(as.numeric(pbiasr:::.AK1logLik(c(0.04, 1, 0.05), d, "log")), expected)
  expect_equal(as.numeric(pbiasr:::.AK2logLik(c(0.04, 1, 1, 1, 0.05), d, "log")), expected)
  expect_equal(pbiasr:::.AK1logLik(c(0.04, 1, 0.05), d, "est"), -sum(expected))
  expect_equal(pbiasr:::.AK2logLik(c(0.04, 1, 1, 1, 0.05), d, "est"), -sum(expected))
  # The ordinary-scale density underflows here, but its log is finite.
  extreme <- d; extreme$effect <- extreme$effect * 1000
  expect_true(is.finite(pbiasr:::.AK1logLik(c(0.04, 1, 0.05), extreme)))
  expect_warning(pbiasr:::.AK1logLik(c(-0.1, 1, 0.05), d), NA)
})

test_that("AK sample tests use the reference degrees of freedom and ML comparison", {
  sample <- pbias_sample_data(1)
  r1 <- pbiasr:::.fit_ak1(sample)
  r2 <- pbiasr:::.fit_ak2(sample)
  expect_equal(r1$status, "ok")
  expect_equal(r2$status, "ok")
  expect_equal(r1$test$df, 1)
  expect_equal(r2$test$df, 3)
  expect_equal(r1$pub_bias_p, 0.4104142, tolerance = 1e-5)
  expect_equal(r2$pub_bias_p, 0.8117845, tolerance = 1e-5)
  for (result in list(r1, r2)) {
    statistic <- 2 * (result$test$loglik_alt - result$test$loglik_null)
    expect_equal(result$pub_bias_p, pchisq(statistic, result$test$df, lower.tail = FALSE))
    b0 <- subset(result$estimates, term == "b0")
    expect_equal(b0$p_value, 2 * pt(-abs(b0$statistic), nrow(sample) - if (result$test$df == 1) 3 else 5))
  }
  # This second sample has positive log-likelihoods, catching the other sign case.
  positive <- pbiasr:::.fit_ak1(pbias_sample_data(2))
  expect_gt(positive$test$loglik_null, 0)
  expect_equal(log(positive$pub_bias_p), log(1.097209e-16), tolerance = 1e-6)
})

test_that("AK tests expose insufficient selection intervals", {
  d <- data.frame(StdID = 1:10, eff = seq(0.01, 0.1, length.out = 10), se = rep(1, 10))
  expect_equal(pbiasr:::.fit_ak1(d)$status, "not_applicable")
  expect_equal(pbiasr:::.fit_ak2(d)$status, "not_applicable")
})
