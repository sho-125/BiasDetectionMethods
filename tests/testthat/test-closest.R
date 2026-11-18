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

test_that("all eligible individual conditions are searched before any summary", {
  x <- pbias_closest(65, .2)
  expect_equal(nrow(x$candidates), 864L)
  expect_equal(x$meta$n_searched, 216L)
  expect_equal(unname(x$meta$searched_counts), c(108L,108L))
  expect_equal(x$matches$Environment, c(292L,484L,309L,485L))
  expect_equal(x$matches$K, c(60,60,100,100))
  expect_equal(x$matches$biased, c(FALSE,TRUE,FALSE,TRUE))
  expect_equal(x$matches$mean_tau, c(.1952589,.2048735,.1968948,.2103092), tolerance=1e-6)
  expect_equal(x$matches$mean_tau2, c(.03935021,.04435287,.03948388,.04563301), tolerance=1e-7)
  expect_true(all(x$matches$n_conditions == 1))
  expect_equal(nrow(x$performance), 84L)
  expect_true(all(is.na(x$performance$balanced_accuracy)))
  expect_equal(x$meta$pooling, "none")
  expect_error(pbias_closest(65,.2,pooling="publication_bias"), "no longer supported")
  data <- readRDS(system.file("extdata","closest_conditions.rds",package="pbiasr"))
  expect_false("profiles" %in% names(data))
  expect_equal(sum(1000 - data$conditions$n_tau), 667L)
  raw <- data$conditions[match(x$matches$Environment,data$conditions$Environment),]
  expect_equal(x$matches$mean_tau,raw$mean_tau)
  expect_equal(x$matches$Censor,raw$Censor)
})

test_that("QRP and Censoring select individual rows independently in both designs", {
  data <- readRDS(system.file("extdata","closest_conditions.rds",package="pbiasr"))
  for (design in c("original","extended")) for (qrp in c("without","with","all")) for (censoring in c("without","with","all")) {
    x <- pbias_closest(65,.2,design=design,qrp=qrp,censoring=censoring)
    expected <- data$conditions
    if(design=="original") expected<-expected[expected$K>=10 & expected$K<=100,]
    if(qrp=="without") expected<-expected[expected$QRP=="none",]
    if(qrp=="with") expected<-expected[expected$QRP!="none",]
    if(censoring=="without") expected<-expected[expected$Censor=="none",]
    if(censoring=="with") expected<-expected[expected$Censor!="none",]
    expect_setequal(x$candidates$Environment,expected$Environment)
    expect_true(all(x$candidates$n_conditions==1))
    for (k in unique(x$matches$K)) {
      possible<-expected$mean_tau[expected$K==k]
      lo<-possible[possible<.2];hi<-possible[possible>.2]
      expect_equal(sort(x$matches$mean_tau[x$matches$K==k]),sort(c(max(lo),min(hi))))
    }
    expect_equal(nrow(x$performance),21*nrow(x$matches))
  }
  pb<-pbias_closest(65,.2,universe="pb_only")
  expect_equal(nrow(pb$candidates),768L)
  expect_true(all(pb$candidates$biased))
  expect_error(pbias_closest(65,.2,qrp="without",censoring="without",universe="pb_only"),"No simulation conditions")
})

test_that("condition rates and post-selection BA match independent raw summaries", {
  x<-pbias_closest(65,.2)
  data<-readRDS(system.file("extdata","closest_conditions.rds",package="pbiasr"))
  for(method in x$summary$method) {
    raw<-data$condition_performance[data$condition_performance$method==method,]
    raw<-raw[match(x$matches$Environment,raw$Environment),]
    actual<-x$performance[x$performance$method==method,]
    rates<-raw$n_rejected/raw$n_valid
    rates[raw$n_valid==0]<-NA_real_
    expect_equal(actual$rejection_rate,rates)
    expect_equal(actual$applicability,raw$n_valid/1000)
    biased<-x$matches$biased
    tpr<-if(all(is.finite(rates[biased])))mean(rates[biased]) else NA_real_
    fpr<-if(all(is.finite(rates[!biased])))mean(rates[!biased]) else NA_real_
    s<-x$summary[x$summary$method==method,]
    expect_equal(s$mean_tpr,tpr)
    expect_equal(s$mean_fpr,fpr)
    expect_equal(s$balanced_accuracy,(tpr+1-fpr)/2)
    expect_equal(c(s$n_biased_matches,s$n_unbiased_matches),c(2L,2L))
  }
  expect_true(all(diff(na.omit(x$summary$balanced_accuracy))<=0))
  expect_output(print(x),"864 eligible conditions; 216 searched")
})

test_that("BA weights classes equally and requires complete rates in both classes", {
  p<-data.frame(method=rep(c("A","B","C"),each=4),match_id=rep(paste0("M",1:4),3),
    rejection_rate=c(.1,.8,.6,1, .2,.4,.6,.8, .1,.8,NA,1),applicability=1)
  s<-pbiasr:::.pbias_closest_summary(p,paste0("M",1:4),c(FALSE,TRUE,TRUE,TRUE))
  expect_equal(s$method,c("A","B","C"))
  expect_equal(s$balanced_accuracy,c(.85,.7,NA))
  expect_equal(s$n_biased_matches,rep(3L,3))
  expect_equal(s$n_unbiased_matches,rep(1L,3))
  expect_equal(s$tpr_available,c(3L,3L,2L))
  only_bias<-pbias_closest(65,.2,qrp="with")
  expect_true(all(is.na(only_bias$summary$balanced_accuracy)))
  expect_true(all(is.na(only_bias$summary$mean_fpr)))
  controls<-pbias_closest(65,.2,qrp="without",censoring="without")
  expect_true(all(is.na(controls$summary$balanced_accuracy)))
  expect_true(all(is.na(controls$summary$mean_tpr)))
})
