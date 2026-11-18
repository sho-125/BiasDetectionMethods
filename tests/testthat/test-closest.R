test_that("bracketing retains exact matches, both sides, ties and boundaries", {
  toy <- expand.grid(K = c(10, 30), mean_tau = c(0.1, 0.3))
  bracket <- pbiasr:::.pbias_bracket
  expect_equal(bracket(toy, 20, 0.2)$mean_tau, c(0.1, 0.3, 0.1, 0.3))
  expect_equal(bracket(toy, 10, 0.2)$K, c(10, 10))
  exact <- bracket(toy, 10, 0.1)
  expect_equal(nrow(exact), 1L)
  expect_identical(exact$heterogeneity_match, "Exact tau")
  expect_equal(nrow(bracket(rbind(toy, toy[1, ]), 10, 0.1)), 2L)
  low <- bracket(toy, 3, 0)
  expect_equal(c(low$K, low$mean_tau), c(10, 0.1))
  high <- bracket(toy, 1000, 1)
  expect_equal(c(high$K, high$mean_tau), c(30, 0.3))
})

test_that("public inputs are validated on the tau scale", {
  for (k in list(NA_real_, Inf, 2, 3.5, numeric(), c(10, 20), "10")) {
    expect_error(pbias_closest(k, 0.2), "k must")
  }
  for (tau in list(NA_real_, Inf, -1, numeric(), c(0.1, 0.2), "0.2")) {
    expect_error(pbias_closest(10, tau), "tau must")
  }
  expect_error(pbias_closest(10, 0.2, design = "invalid"), "arg")
})

test_that("individual conditions remain unpooled and rates retain their denominators", {
  data <- readRDS(system.file("extdata", "closest_conditions.rds", package = "pbiasr"))
  expect_equal(nrow(data$conditions), 864L)
  expect_equal(sum(1000 - data$conditions$n_tau), 667L)
  expect_true(all(data$conditions$n_tau > 0))
  expect_equal(nrow(data$condition_performance), 864L * 21L)
  for (design in c("original", "extended")) {
    for (qrp in c("without", "with", "all")) {
      for (universe in c("pb_only", "all")) {
        out <- pbias_closest(65, 0.2, design, qrp, universe)
        base_count <- if (design == "original") 432L else 864L
        expected <- base_count * switch(qrp, without = 1/3, with = 2/3, all = 1) *
          if (universe == "pb_only") 2/3 else 1
        expect_equal(nrow(out$candidates), expected)
        expect_true(all(out$matches$n_conditions == 1))
        expect_equal(sort(unique(out$matches$K)), c(60, 100))
        expect_equal(nrow(out$performance), nrow(out$matches) * 21L)
        expect_true(all(is.na(out$performance$balanced_accuracy)))
        raw <- data$conditions[match(out$matches$Environment, data$conditions$Environment), ]
        expect_equal(out$matches$mean_tau, raw$mean_tau)
        expect_equal(out$matches$Censor, raw$Censor)
        perf <- out$performance
        lookup <- match(paste(out$matches$Environment[match(perf$match_id, out$matches$match_id)], perf$method),
                        paste(data$condition_performance$Environment, data$condition_performance$method))
        ref <- data$condition_performance[lookup, ]
        expect_equal(perf$rejection_rate, ref$rate)
        expect_equal(perf$applicability, ref$n_valid / 1000)
        expect_equal(perf$n_total, rep(1000L, nrow(perf)))
        if (universe == "pb_only") expect_true(all(is.na(perf$fpr)))
      }
    }
  }
  all_conditions <- pbias_closest(65, 0.2, qrp = "all", universe = "all")$candidates
  ref <- subset(all_conditions, K == 30 & true_effect == 0.2 & population_tau == 0.2 & QRP == "none")
  exact <- pbias_closest(30, ref$mean_tau[ref$Censor == "none"], universe = "all")
  expect_true(all(is.na(exact$performance$tpr)))
  expect_equal(exact$performance$fpr, exact$performance$rejection_rate)
})

test_that("optional pooling is explicit and matches independent component averages", {
  pb <- pbias_closest(65, 0.2, qrp = "all", pooling = "publication_bias")
  all <- pbias_closest(65, 0.2, qrp = "all", universe = "all", pooling = "publication_bias")
  expect_setequal(pb$candidates$base_id, all$candidates$base_id)
  expect_true(all(pb$candidates$n_conditions == 2L))
  expect_true(all(all$candidates$n_conditions == 3L))
  data <- readRDS(system.file("extdata", "closest_conditions.rds", package = "pbiasr"))
  for (i in seq_len(nrow(all$matches))) {
    m <- all$matches[i, ]
    parts <- data$conditions[data$conditions$base_id == m$base_id, ]
    expect_equal(m$mean_tau, mean(parts$mean_tau), tolerance = 1e-14)
    expect_gt(abs(m$mean_tau - sqrt(mean(parts$mean_tau2))), 1e-6)
    raw <- subset(data$condition_performance, Environment %in% parts$Environment & method == "TESS")
    result <- subset(all$performance, match_id == m$match_id & method == "TESS")
    expect_equal(result$rejection_rate, mean(raw$rate), tolerance = 1e-14)
    expect_equal(result$n_valid, sum(raw$n_valid))
  }
  expect_true(all(is.na(pb$performance$fpr)))
  expect_output(print(pb), "Optional pooling")
  expect_output(print(pbias_closest(65, 0.2)), "Individual conditions")
  expect_length(pbias_closest(1000, 10)$notes, 2L)
})
