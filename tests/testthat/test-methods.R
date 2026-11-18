test_that("constant standard errors cannot masquerade as Egger slope tests", {
  d <- data.frame(StdID = 1:12, eff = seq(-0.4, 0.6, length.out = 12), se = rep(0.2, 12))
  for (fit in list(pbiasr:::.fit_egger_fe, pbiasr:::.fit_egger_re, pbiasr:::.fit_egger_wls)) {
    result <- fit(d)
    expect_equal(result$status, "not_applicable")
    expect_true(is.na(result$pub_bias_p))
  }
  result <- suppressWarnings(pbias_table(d$eff, d$se))
  rows <- result$bias_summary[result$bias_summary$paper_method %in% c("EGGER1", "EGGER2", "EGGER3"), ]
  expect_true(all(rows$status == "not_applicable"))
  expect_true(all(is.na(rows[["p-value"]])))
})

test_that("unidentified sample-size instruments are marked unavailable", {
  d <- data.frame(StdID = 1:12, eff = seq(-0.4, 0.6, length.out = 12),
                  se = seq(0.1, 0.3, length.out = 12), obs = rep(100, 12))
  expect_equal(pbiasr:::.fit_maive(d)$status, "not_applicable")
  expect_equal(pbiasr:::.fit_fativ(d)$status, "not_applicable")
})

test_that("sample runs return exactly 21 bounded p-values or explicit unavailable states", {
  for (i in 1:2) {
    d <- pbias_sample_data(i)
    result <- suppressWarnings(pbias_table(d$eff, d$se, obs = d$obs))
    b <- result$bias_summary
    expect_equal(nrow(b), 21L)
    expect_equal(length(unique(b$paper_method)), 21L)
    p <- b[["p-value"]]
    expect_true(all(is.na(p) | (is.finite(p) & p >= 0 & p <= 1)))
    expect_false(any(b$status == "ok" & is.na(p)))
    expect_equal(nrow(result$publication_table), 21L)
  }
  expect_error(pbias_table(1:3, rep(0.2, 3), alpha = NA_real_), "single number")
})

test_that("bundled source uses its own sample data and performance resources", {
  expect_true(file.exists(pbiasr:::.pbias_extdata_file("SampleData1.RData")))
  performance <- pbiasr:::.type1_performance_table("SS")
  expect_gt(nrow(performance), 0L)
  expect_true(all(performance$method %in% pbiasr:::.method_map()$paper_method))
  expect_named(pbias_sample_data("all"), c("SampleData1", "SampleData2"))
})

test_that("TES and caliper thresholds reproduce the simulation implementation", {
  sample1 <- pbias_sample_data(1)
  sample2 <- pbias_sample_data(2)
  expect_equal(pbiasr:::.fit_tes(sample1)$pub_bias_p, 0.71287493, tolerance = 1e-6)
  expect_equal(log(pbiasr:::.fit_tes(sample2)$pub_bias_p), log(1.315190e-9), tolerance = 1e-6)
  result <- suppressWarnings(pbias_table(sample1$eff, sample1$se, obs = sample1$obs))
  tes <- subset(result$bias_summary, paper_method == "TES")
  expect_equal(tes[["p-value"]], 0.71287493, tolerance = 1e-6)
  calipers <- pbiasr:::.fit_calipers(sample1)$summary
  expect_equal(calipers$publication_bias_p[3:4], c(0.9375, 0.91021728515625))
  # These values lie in the rounded CALI05 window but outside +/-5% of qnorm(.975).
  boundary <- data.frame(eff = c(rep(1.861, 3), rep(2.059, 3), 1.86, 2.06, 1.96), se = 1)
  result <- pbiasr:::.fit_calipers(boundary)$summary
  expect_equal(result$note[1], "Counts: over=3, under=3")
  expect_equal(result$publication_bias_p[1], binom.test(3, 6, alternative = "greater")$p.value)
})
