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
        out <- pbias_closest(65, 0.2, design, qrp, universe, pooling = "none")
        base_count <- if (design == "original") 432L else 864L
        expected <- base_count * switch(qrp, without = 1/3, with = 2/3, all = 1) *
          if (universe == "pb_only") switch(qrp, without = 2/3, with = 1, all = 8/9) else 1
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
  all_conditions <- pbias_closest(65, 0.2, qrp = "all", universe = "all", pooling = "none")$candidates
  ref <- subset(all_conditions, K == 30 & true_effect == 0.2 & population_tau == 0.2 & QRP == "none")
  exact <- pbias_closest(30, ref$mean_tau[ref$Censor == "none"], universe = "all", qrp = "without", pooling = "none")
  expect_true(all(is.na(exact$performance$tpr)))
  expect_equal(exact$performance$fpr, exact$performance$rejection_rate)
})

test_that("optional pooling is explicit and matches independent component averages", {
  pb <- pbias_closest(65, 0.2, qrp = "all", universe = "pb_only", pooling = "publication_bias")
  all <- pbias_closest(65, 0.2, qrp = "all", universe = "all", pooling = "publication_bias")
  expect_setequal(pb$candidates$base_id, all$candidates$base_id)
  expect_equal(pb$candidates$n_conditions, ifelse(pb$candidates$QRP == "none", 2L, 3L))
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
  expect_output(print(pb), "Censoring: all")
  expect_output(print(pbias_closest(65, 0.2)), "Mean_tau2")
  expect_length(pbias_closest(1000, 10)$notes, 2L)
})


test_that("new defaults are grouped Extended / all QRP / All and retain both tau scales", {
  x <- pbias_closest(65, 0.2)
  expect_equal(x$meta[c("design", "qrp", "universe", "pooling")],
               list(design = "extended", qrp = "all", universe = "all", pooling = "publication_bias"))
  expect_equal(nrow(x$candidates), 288L)
  expect_equal(nrow(x$matches), 4L)
  expect_equal(x$matches$K, c(60, 60, 100, 100))
  expect_equal(x$matches$mean_tau, c(.1811738, .2175545, .1842785, .2197795), tolerance = 1e-6)
  expect_equal(x$matches$mean_tau2, c(.03393001, .04971837, .03462290, .04996398), tolerance = 1e-7)
  expect_true(all(abs(x$matches$mean_tau2 - x$matches$mean_tau^2) > 1e-5))
  expect_equal(nrow(x$summary), 21L)
  expect_true(all(diff(na.omit(x$summary$mean_balanced_accuracy)) <= 0))
  for (method in x$summary$method) {
    b <- x$performance[x$performance$method == method, ]$balanced_accuracy
    expected <- if (all(is.finite(b))) mean(b) else NA_real_
    expect_equal(x$summary$mean_balanced_accuracy[x$summary$method == method], expected)
  }
})

test_that("biased-only membership includes QRP-only conditions in raw and pooled modes", {
  all <- pbias_closest(65, .2, pooling = "none")
  pb <- pbias_closest(65, .2, universe = "pb_only", pooling = "none")
  expected <- subset(all$candidates, Censor != "none" | QRP != "none")
  expect_setequal(pb$candidates$Environment, expected$Environment)
  expect_equal(nrow(pb$candidates), 768L)
  expect_equal(sum(pb$candidates$Censor == "none"), 192L)
  with_qrp <- pbias_closest(65, .2, qrp = "with", universe = "pb_only", pooling = "none")
  expect_equal(nrow(with_qrp$candidates), 576L)
  pooled_pb <- pbias_closest(65, .2, universe = "pb_only")$candidates
  pooled_all <- pbias_closest(65, .2)$candidates
  positive <- pooled_pb$QRP != "none"
  i <- match(pooled_pb$base_id[positive], pooled_all$base_id)
  expect_equal(pooled_pb$mean_tau[positive], pooled_all$mean_tau[i])
  expect_equal(pooled_pb$mean_tau2[positive], pooled_all$mean_tau2[i])
  expect_true(all(pooled_pb$n_conditions[positive] == 3))
})

test_that("summary gives equal match weights and never fills missing BA", {
  p <- data.frame(method = rep(c("A", "B", "C"), each = 2), match_id = rep(c("M1", "M2"), 3),
                  balanced_accuracy = c(.6, .8, .4, .6, .7, NA),
                  rejection_rate = c(.2, .4, .8, .6, .3, .5), applicability = c(1, 1, .8, 1, 1, 1))
  s <- pbiasr:::.pbias_closest_summary(p, c("M1", "M2"))
  expect_equal(s$method, c("A", "B", "C"))
  expect_equal(s$mean_balanced_accuracy, c(.7, .5, NA))
  expect_equal(s$mean_rejection_rate, c(.3, .7, .4))
  expect_equal(s$ba_available, c(2L, 2L, 1L))
  expect_equal(s$low_applicability_matches, c(0L, 1L, 0L))
})


test_that("independent QRP and Censoring filters use exactly the requested components", {
  data <- readRDS(system.file("extdata", "closest_conditions.rds", package = "pbiasr"))
  for (design in c("original", "extended")) {
    for (qrp in c("without", "with", "all")) {
      for (censoring in c("without", "with", "all")) {
        x <- pbias_closest(65, .2, design = design, qrp = qrp, censoring = censoring)
        expect_equal(x$meta$censoring, censoring)
        expect_equal(unique(x$candidates$n_conditions), switch(censoring, without = 1L, with = 2L, all = 3L))
        expect_equal(nrow(x$candidates), (if (design == "original") 144L else 288L) * switch(qrp, without = 1/3, with = 2/3, all = 1))
        for (i in seq_len(nrow(x$matches))) {
          m <- x$matches[i, ]
          parts <- data$conditions[data$conditions$base_id == m$base_id, ]
          if (censoring == "without") parts <- parts[parts$Censor == "none", ]
          if (censoring == "with") parts <- parts[parts$Censor != "none", ]
          expect_equal(m$mean_tau, mean(parts$mean_tau), tolerance = 1e-14)
          expect_equal(m$mean_tau2, mean(parts$mean_tau2), tolerance = 1e-14)
          expect_equal(m$n_biased_conditions, sum(parts$biased))
          expect_equal(m$n_unbiased_conditions, sum(!parts$biased))
          raw <- subset(data$condition_performance, Environment %in% parts$Environment & method == "BEGG")
          actual <- subset(x$performance, match_id == m$match_id & method == "BEGG")
          expect_equal(actual$rejection_rate, mean(raw$rate), tolerance = 1e-14)
        }
      }
    }
  }
  qrp_only <- pbias_closest(65, .2, qrp = "with", censoring = "without")
  expect_true(all(qrp_only$candidates$n_biased_conditions == 1))
  expect_true(all(qrp_only$candidates$n_unbiased_conditions == 0))
  expect_true(all(is.na(qrp_only$performance$balanced_accuracy)))
  control <- pbias_closest(65, .2, qrp = "without", censoring = "without")
  expect_true(all(control$candidates$n_unbiased_conditions == 1))
  expect_error(pbias_closest(65, .2, censoring = "invalid"), "arg")
  expect_error(pbias_closest(65, .2, qrp = "without", censoring = "without", universe = "pb_only"), "No simulation conditions")
})
