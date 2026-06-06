
# User-facing function -----------------------------------------------------

#' Run publication-bias methods for a meta-analysis
#'
#' @param yi Numeric vector of observed effect sizes.
#' @param sei Numeric vector of standard errors associated with \code{yi}.
#' @param obs Optional numeric vector of study-level sample size proxies. This is
#'   required for the EGGER4 / EGGER5 implementations in this package.
#' @param study_id Optional study identifiers. Defaults to 1:n.
#' @param alpha Significance level used by the caliper windows and excess-significance
#'   tests. The core implementations follow the uploaded code and paper, which use
#'   0.05 by default.
#' @param sort_by_rank Logical. If \code{TRUE}, sort the displayed publication-bias
#'   results by the Balanced Accuracy ranking for the detected paper case.
#'
#' @return An object of class \code{pbias_result} with elements:
#'   \itemize{
#'     \item \code{bias_summary}: one row per paper method, with publication-bias p-values.
#'     \item \code{publication_table}: formatted display table of publication-bias tests.
#'     \item \code{case}: the paper case based on number of studies and I-squared.
#'     \item \code{interpretation}: short interpretation text for the publication-bias tests.
#'     \item \code{meta}: run metadata.
#'   }
#'
#' @examples
#' yi <- c(0.20, 0.35, 0.10, 0.42, 0.30, 0.28)
#' sei <- c(0.10, 0.12, 0.15, 0.11, 0.09, 0.13)
#' out <- pbias_table(yi, sei)
#' print(out)
#'
#' @export
pbias_table <- function(yi, sei, obs = NULL, study_id = NULL, alpha = 0.05,
                        sort_by_rank = TRUE) {
  .validate_inputs(yi, sei, obs, study_id, alpha, sort_by_rank)

  n <- length(yi)
  if (is.null(study_id)) {
    study_id <- seq_len(n)
  }

  dt <- data.frame(
    eff = as.numeric(yi),
    se = as.numeric(sei),
    StdID = as.numeric(study_id),
    stringsAsFactors = FALSE
  )
  if (!is.null(obs)) {
    dt$obs <- as.numeric(obs)
  }

  # Metadata and mappings
  method_map <- .method_map()
  bias_rows <- vector("list", nrow(method_map))
  names(bias_rows) <- method_map$paper_method

  # Helper for summary rows
  add_bias <- function(code_method, p_value, status = "ok", note = "") {
    idx <- match(code_method, method_map$code_method)
    meta_row <- method_map[idx, , drop = FALSE]
    bias_rows[[meta_row$paper_method]] <<- data.frame(
      paper_method = meta_row$paper_method,
      code_method = meta_row$code_method,
      category = meta_row$category,
      publication_bias_p = as.numeric(p_value),
      status = status,
      note = note,
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  # Baseline random-effects model (not counted among the 21 paper methods)
  baseline <- .fit_baseline_re(dt)
  case_info <- .classify_meta_case(n, baseline$model)
  performance_rank <- .type1_performance_table(case_info$code)

  # EGGER1 - FE
  res <- .fit_egger_fe(dt)
  add_bias("FE", res$pub_bias_p, res$status, res$note)

  # EGGER2 - RE
  res <- .fit_egger_re(dt)
  add_bias("RE", res$pub_bias_p, res$status, res$note)

  # EGGER3 - WLS
  res <- .fit_egger_wls(dt)
  add_bias("WLS", res$pub_bias_p, res$status, res$note)

  # EGGER4 - MAIVE
  res <- .fit_maive(dt)
  add_bias("MAIVE", res$pub_bias_p, res$status, res$note)

  # EGGER5 - FATIV
  res <- .fit_fativ(dt)
  add_bias("FATIV", res$pub_bias_p, res$status, res$note)

  # Selection models
  res <- .fit_psm3(dt)
  add_bias("3PSM", res$pub_bias_p, res$status, res$note)

  res <- .fit_psm4(dt)
  add_bias("4PSM", res$pub_bias_p, res$status, res$note)

  res <- .fit_ak1(dt)
  add_bias("AK1", res$pub_bias_p, res$status, res$note)

  res <- .fit_ak2(dt)
  add_bias("AK2", res$pub_bias_p, res$status, res$note)

  res <- .fit_puniform(dt)
  add_bias("Puniform", res$pub_bias_p, res$status, res$note)

  # Other SE-based methods
  res <- .fit_ek(dt)
  add_bias("EK", res$pub_bias_p, res$status, res$note)

  res <- .fit_begg(dt)
  add_bias("Begg", res$pub_bias_p, res$status, res$note)

  res <- .fit_skew(dt, baseline$model)
  add_bias("SKEWNESS", res$skew_p, res$status, res$note)
  add_bias("SKEWCombined", res$combined_p, res$status, res$note)

  # Excess significance tests
  res <- .fit_tes(dt, alpha = alpha)
  add_bias("TES", res$pub_bias_p, res$status, res$note)

  res <- .fit_psst_tess(dt)
  add_bias("PSST", res$psst_p, res$status, res$note)
  add_bias("TESS", res$tess_p, res$status, res$note)

  # Caliper tests
  res <- .fit_calipers(dt, alpha = alpha)
  for (i in seq_len(nrow(res$summary))) {
    add_bias(res$summary$code_method[i], res$summary$publication_bias_p[i],
             res$summary$status[i], res$summary$note[i])
  }

  # Fill any missing summary rows with NA
  for (i in seq_len(nrow(method_map))) {
    if (is.null(bias_rows[[method_map$paper_method[i]]])) {
      add_bias(method_map$code_method[i], NA_real_, "not_run", "No result returned.")
    }
  }

  bias_summary <- do.call(rbind, bias_rows)
  row.names(bias_summary) <- NULL

  out <- list(
    bias_summary = .make_bias_summary(bias_summary),
    publication_table = .make_publication_table(
      bias_summary,
      performance_rank = performance_rank,
      sort_by_rank = sort_by_rank
    ),
    case = case_info,
    interpretation = paste(
      "Null hypothesis: no publication bias.",
      "Rejection indicates that the method detects the existence of publication bias."
    ),
    meta = list(
      n_studies = n,
      has_obs = !is.null(obs),
      alpha = alpha,
      sort_by_rank = sort_by_rank,
      ranking_metric = "Balanced Accuracy",
      performance_source = "Performance tables from the publication bias paper",
      timestamp = Sys.time()
    )
  )
  class(out) <- "pbias_result"
  out
}

#' @export
print.pbias_result <- function(x, digits = 3, ...) {
  cat("\n")
  cat("Bias Test Results\n")
  cat("=============================\n")
  cat("Caution: These recommendations are derived from simulated meta-analyses based on the experimental conditions in Carter et al. (2019). Rankings reflect Balanced Accuracy in those simulations, not universal method performance, and should be applied with care because bias-detection tests may behave differently in real meta-analytic applications.\n\n")
  cat("Null hypothesis: no bias\n")
  cat("\n")
  cat("Studies       : ", x$meta$n_studies, "\n", sep = "")
  cat("obs supplied  : ", if (isTRUE(x$meta$has_obs)) "yes" else "no", "\n", sep = "")
  if (!is.null(x$case)) {
    cat("Paper case    : ", .format_case_value(x$case, digits = digits), "\n", sep = "")
  }
  cat("\n")
  cat("Ranked by Balanced Accuracy\n")
  cat("---------------------------\n")
  .print_publication_table(x$publication_table)
  cat("\n")
  cat("-----\n")
  cat(paste(strwrap(.format_output_note(x), width = 78), collapse = "\n"), "\n")
  cat("\n")

  invisible(x)
}

# Internal maps and templates ---------------------------------------------

.method_map <- function() {
  data.frame(
    paper_method = c(
      "BEGG", "EGGER1", "EGGER2", "EGGER3", "EGGER4", "EGGER5",
      "EK", "SKEW1", "SKEW2", "PSM3", "PSM4", "AK1", "AK2", "PUNIF",
      "TES", "PSST", "TESS", "CALI05", "CALI10", "CALI15", "CALI20"
    ),
    code_method = c(
      "Begg", "FE", "RE", "WLS", "MAIVE", "FATIV",
      "EK", "SKEWNESS", "SKEWCombined", "3PSM", "4PSM", "AK1", "AK2", "Puniform",
      "TES", "PSST", "TESS", "Caliper05", "Caliper10", "Caliper15", "Caliper20"
    ),
    category = c(
      rep("standard_error_based", 9),
      rep("selection_model", 5),
      rep("excess_statistical_significance", 3),
      rep("caliper", 4)
    ),
    stringsAsFactors = FALSE
  )
}

.validate_inputs <- function(yi, sei, obs, study_id, alpha, sort_by_rank) {
  if (!is.numeric(yi) || !is.numeric(sei)) {
    stop("`yi` and `sei` must be numeric vectors.")
  }
  if (length(yi) != length(sei)) {
    stop("`yi` and `sei` must have the same length.")
  }
  if (length(yi) < 3L) {
    stop("At least 3 studies are required.")
  }
  if (any(!is.finite(yi)) || any(!is.finite(sei)) || any(sei <= 0)) {
    stop("`yi` must be finite and `sei` must be finite and strictly positive.")
  }
  if (!is.null(obs)) {
    if (!is.numeric(obs) || length(obs) != length(yi)) {
      stop("If supplied, `obs` must be a numeric vector of the same length as `yi`.")
    }
    if (any(!is.finite(obs) | obs <= 0)) {
      stop("If supplied, `obs` must be finite and strictly positive.")
    }
  }
  if (!is.null(study_id) && length(study_id) != length(yi)) {
    stop("If supplied, `study_id` must have the same length as `yi`.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number between 0 and 1.")
  }
  if (!is.logical(sort_by_rank) || length(sort_by_rank) != 1L || is.na(sort_by_rank)) {
    stop("`sort_by_rank` must be TRUE or FALSE.")
  }
}

.na_result <- function(note = "", status = "failed") {
  list(
    estimates = NULL,
    pub_bias_p = NA_real_,
    status = status,
    note = note
  )
}

.make_estimates <- function(term, estimate = NA_real_, std_error = NA_real_,
                            statistic = NA_real_, p_value = NA_real_,
                            conf_low = NA_real_, conf_high = NA_real_,
                            note = "") {
  data.frame(
    term = term,
    estimate = as.numeric(estimate),
    std_error = as.numeric(std_error),
    statistic = as.numeric(statistic),
    p_value = as.numeric(p_value),
    conf_low = as.numeric(conf_low),
    conf_high = as.numeric(conf_high),
    note = as.character(note),
    stringsAsFactors = FALSE
  )
}

.format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

.format_result_p <- function(p, status, digits = 3) {
  if (identical(status, "ok") && !is.na(p)) return(.format_num(p, digits))
  if (identical(status, "requires_obs")) return("requires obs")
  if (identical(status, "not_applicable")) return("n/a")
  if (identical(status, "failed")) return("failed")
  if (identical(status, "not_run")) return("not run")
  if (!is.na(p)) return(.format_num(p, digits))
  ""
}

.format_output_note <- function(x) {
  case_code <- x$case$code
  
  case_text <- if (!is.null(case_code) && length(case_code) == 1 && !is.na(case_code)) {
    paste0("case ", case_code)
  } else {
    "the detected paper case"
  }
  
  paste(
    "* Some tests were excluded from the rankings because the data conditions did not satisfy their applicability requirements or because they exhibited substantial convergence problems.\n",
    "Note: The estimated p-value is calculated from the user's data.",
    paste0(
      "Balanced Accuracy and Logit Distance are based on the results from OUR BIAS paper for ",
      case_text,
      "."
    ),
    "Rows are sorted by Balanced Accuracy when sort_by_rank = TRUE.",
    "A small p-value indicates evidence of bias.",
    sep = "\n"
  )
}

.print_publication_table <- function(tab) {
  if (is.null(tab) || nrow(tab) == 0) {
    cat("(no publication-bias results)\n")
    return(invisible(NULL))
  }

  work <- as.data.frame(tab, stringsAsFactors = FALSE, check.names = FALSE)
  for (nm in names(work)) {
    work[[nm]] <- as.character(work[[nm]])
    work[[nm]][is.na(work[[nm]])] <- ""
  }

  widths <- vapply(seq_along(work), function(i) {
    max(nchar(c(names(work)[i], work[[i]]), type = "width"))
  }, integer(1))

  header <- paste(mapply(function(value, width) {
    format(value, width = width, justify = "left")
  }, names(work), widths, USE.NAMES = FALSE), collapse = "  ")
  rule <- paste(vapply(widths, function(w) paste(rep("-", w), collapse = ""),
                       character(1)), collapse = "  ")
  cat(header, "\n", sep = "")
  cat(rule, "\n", sep = "")

  for (i in seq_len(nrow(work))) {
    row_values <- unlist(work[i, ], use.names = FALSE)
    row <- paste(mapply(function(value, width) {
      format(value, width = width, justify = "left")
    }, row_values, widths, USE.NAMES = FALSE), collapse = "  ")
    cat(row, "\n", sep = "")
  }

  invisible(NULL)
}

.classify_meta_case <- function(n, model) {
  i2 <- NA_real_
  if (!is.null(model) &&
      !is.null(model$I2) &&
      length(model$I2) == 1L &&
      is.finite(model$I2)) {
    i2 <- max(0, as.numeric(model$I2) / 100)
  }

  if (n <= 10L) {
    size_code <- "S"
    size_category <- "Small"
    size_rule <- "K <= 10"
  } else if (n <= 100L) {
    size_code <- "M"
    size_category <- "Medium"
    size_rule <- "10 < K <= 100"
  } else {
    size_code <- "L"
    size_category <- "Large"
    size_rule <- "K > 100"
  }

  if (!is.finite(i2)) {
    heterogeneity_code <- NA_character_
    heterogeneity_category <- "Unknown"
    heterogeneity_rule <- "I2 unavailable"
    code <- NA_character_
    note <- "Random-effects model failed; I2-based paper case could not be assigned."
  } else if (i2 < 0.29) {
    heterogeneity_code <- "S"
    heterogeneity_category <- "Small"
    heterogeneity_rule <- "I2 < 0.29"
    code <- paste0(size_code, heterogeneity_code)
    note <- ""
  } else if (i2 < 0.60) {
    heterogeneity_code <- "M"
    heterogeneity_category <- "Medium"
    heterogeneity_rule <- "0.29 <= I2 < 0.60"
    code <- paste0(size_code, heterogeneity_code)
    note <- ""
  } else {
    heterogeneity_code <- "L"
    heterogeneity_category <- "Large"
    heterogeneity_rule <- "I2 >= 0.60"
    code <- paste0(size_code, heterogeneity_code)
    note <- ""
  }

  list(
    code = code,
    n_studies = n,
    i2 = i2,
    size_code = size_code,
    size_category = size_category,
    size_rule = size_rule,
    heterogeneity_code = heterogeneity_code,
    heterogeneity_category = heterogeneity_category,
    heterogeneity_rule = heterogeneity_rule,
    note = note
  )
}

.format_case_value <- function(case, digits = 3) {
  if (is.null(case)) {
    return("unavailable")
  }
  if (is.na(case$code)) {
    return(paste0("unavailable; K=", case$n_studies, "; I2 unavailable"))
  }
  paste0(
    case$code,
    " (", case$size_category, " size: ", case$size_rule,
    "; ", case$heterogeneity_category, " heterogeneity: ", case$heterogeneity_rule,
    "; I2=", .format_num(case$i2, digits), ")"
  )
}

.type1_performance_table <- function(case_code = NULL) {
  if (is.null(case_code) || is.na(case_code)) {
    return(data.frame())
  }

  path <- system.file("extdata", "type1_performance.csv", package = "pbiasr")
  if (!nzchar(path)) {
    local_path <- file.path("inst", "extdata", "type1_performance.csv")
    if (file.exists(local_path)) {
      path <- local_path
    }
  }
  if (!nzchar(path) || !file.exists(path)) {
    return(.type1_performance_fallback(case_code))
  }

  tab <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  num_cols <- setdiff(names(tab), c("case", "method"))
  for (nm in num_cols) {
    tab[[nm]] <- as.numeric(tab[[nm]])
  }
  tab <- tab[tab$case == case_code, , drop = FALSE]
  if (nrow(tab) == 0) {
    return(.type1_performance_fallback(case_code))
  }
  .rank_performance_table(tab)
}

.rank_performance_table <- function(tab) {
  if (is.null(tab) || nrow(tab) == 0) {
    return(data.frame())
  }
  tab$rank <- as.integer(rank(-tab$balanced_accuracy, ties.method = "min"))
  tab <- tab[order(tab$rank, tab$false_positive, tab$method), , drop = FALSE]
  tab[, c("case", "rank", "method", "true_positive", "false_positive",
          "precision", "sensitivity", "f1", "balanced_accuracy", "inv_fpr", "ldist")]
}

.make_publication_table <- function(bias_summary, performance_rank = NULL,
                                    sort_by_rank = TRUE, digits = 3) {
  pvalue <- mapply(
    .format_result_p,
    bias_summary$publication_bias_p,
    bias_summary$status,
    MoreArgs = list(digits = digits),
    USE.NAMES = FALSE
  )
  
  tab <- data.frame(
    method = bias_summary$paper_method,
    pvalue = pvalue,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!is.null(performance_rank) && nrow(performance_rank) > 0) {
    idx <- match(tab$method, performance_rank$method)
    
    rank <- performance_rank$rank[idx]
    
    balanced_accuracy <- performance_rank$balanced_accuracy[idx]
    balanced_accuracy_text <- ifelse(
      is.na(balanced_accuracy),
      "Unranked*",
      .format_num(balanced_accuracy, digits)
    )
    
    ldist <- performance_rank$ldist[idx]
    ldist_text <- ifelse(
      is.na(ldist),
      "",
      .format_num(ldist, digits)
    )
    
    tab <- data.frame(
      Method = tab$method,
      `P-value` = tab$pvalue,
      `Balanced Accuracy` = balanced_accuracy_text,
      `Logit Distance`= ldist_text,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    if (isTRUE(sort_by_rank)) {
      order_key <- ifelse(is.na(rank), Inf, rank)
      tab <- tab[order(order_key, seq_len(nrow(tab))), , drop = FALSE]
    }
  } else {
    tab <- data.frame(
      Method = tab$method,
      `P-value` = tab$pvalue,
      `Balanced Accuracy` = "not available",
      `Logit Distance` = "not available",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  tab
}

.make_bias_summary <- function(bias_summary) {
  out <- bias_summary[, setdiff(names(bias_summary), "code_method"), drop = FALSE]
  names(out)[names(out) == "publication_bias_p"] <- "p-value"
  row.names(out) <- NULL
  out
}

# Baseline ----------------------------------------------------------------

.fit_baseline_re <- function(dt) {
  mod <- tryCatch(
    metafor::rma(dt$eff, dt$se^2, method = "ML"),
    error = function(e) NULL
  )
  if (is.null(mod)) {
    return(list(
      model = NULL,
      estimates = rbind(
        .make_estimates("b0", note = "Baseline RE model failed."),
        .make_estimates("tau2"),
        .make_estimates("I2")
      )
    ))
  }

  est <- rbind(
    .make_estimates("b0", mod$b[[1]], mod$se[[1]], mod$zval[[1]], mod$pval[[1]],
                    mod$ci.lb[[1]], mod$ci.ub[[1]]),
    .make_estimates("tau2", mod$tau2),
    .make_estimates("I2", mod$I2 / 100)
  )
  list(model = mod, estimates = est)
}

# Egger variants -----------------------------------------------------------

.fit_egger_fe <- function(dt) {
  mod <- tryCatch(
    metafor::rma(dt$eff ~ dt$se, dt$se^2, method = "FE"),
    error = function(e) NULL
  )
  if (is.null(mod)) return(.na_result("FE Egger model failed."))

  est <- .make_estimates(
    term = c("b0", "se"),
    estimate = as.numeric(mod$b),
    std_error = mod$se,
    statistic = mod$zval,
    p_value = mod$pval,
    conf_low = mod$ci.lb,
    conf_high = mod$ci.ub
  )
  list(estimates = est, pub_bias_p = est$p_value[est$term == "se"], status = "ok", note = "")
}

.fit_egger_re <- function(dt) {
  mod <- tryCatch(
    metafor::rma(dt$eff ~ dt$se, dt$se^2, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(mod)) return(.na_result("RE Egger model failed."))

  est <- .make_estimates(
    term = c("b0", "se"),
    estimate = as.numeric(mod$b),
    std_error = mod$se,
    statistic = mod$zval,
    p_value = mod$pval,
    conf_low = mod$ci.lb,
    conf_high = mod$ci.ub
  )
  est <- rbind(est, .make_estimates("tau2", mod$tau2))
  list(estimates = est, pub_bias_p = est$p_value[est$term == "se"][1], status = "ok", note = "")
}

.fit_egger_wls <- function(dt) {
  mod <- tryCatch(
    stats::lm(dt$eff ~ dt$se, weights = 1 / (dt$se^2)),
    error = function(e) NULL
  )
  if (is.null(mod)) return(.na_result("WLS Egger model failed."))

  co <- summary(mod)$coefficients
  ci <- tryCatch(stats::confint(mod), error = function(e) matrix(NA_real_, nrow = 2, ncol = 2))
  est <- .make_estimates(
    term = c("b0", "se"),
    estimate = as.numeric(co[, 1]),
    std_error = as.numeric(co[, 2]),
    statistic = as.numeric(co[, 3]),
    p_value = as.numeric(co[, 4]),
    conf_low = as.numeric(ci[, 1]),
    conf_high = as.numeric(ci[, 2])
  )
  list(estimates = est, pub_bias_p = est$p_value[est$term == "se"], status = "ok", note = "")
}

# MAIVE / FATIV ------------------------------------------------------------

.MAIVE_EST <- function(dt) {
  invObs <- 1 / dt$obs
  se2 <- dt$se^2
  fit.se2 <- stats::lm(se2 ~ 1 + invObs)
  if (fit.se2$coefficients[1] < 0) {
    fit.se2 <- stats::lm(se2 ~ 0 + invObs)
  }
  se2fit <- stats::fitted(fit.se2)
  se1fit <- sqrt(se2fit)
  maive.est <- stats::lm(dt$eff ~ 1 + se1fit)
  fatpet.test <- summary(maive.est)$coefficient[2, 4]
  if (summary(maive.est)$coefficient[1, 4] < 0.05) {
    maive.est <- stats::lm(dt$eff ~ 1 + se2fit)
  }
  list(model = maive.est, fatpet_test = fatpet.test)
}

.MAIVE2_EST <- function(dt) {
  invObs <- 1 / dt$obs
  se2 <- dt$se^2
  fit.se2 <- stats::lm(se2 ~ 1 + invObs)
  if (fit.se2$coefficients[1] < 0) {
    fit.se2 <- stats::lm(se2 ~ 0 + invObs)
  }
  se2fit <- stats::fitted(fit.se2)
  se1fit <- sqrt(se2fit)
  maive.est <- stats::lm(dt$eff ~ 1 + se1fit, weights = 1 / se2)
  fatpet.test <- summary(maive.est)$coefficient[2, 4]
  if (summary(maive.est)$coefficient[1, 4] < 0.05) {
    maive.est <- stats::lm(dt$eff ~ 1 + se2fit, weights = 1 / se2)
  }
  list(model = maive.est, fatpet_test = fatpet.test)
}

.fit_maive <- function(dt) {
  if (is.null(dt$obs)) {
    return(.na_result("EGGER4 / MAIVE requires `obs`.", status = "requires_obs"))
  }
  ans <- tryCatch(.MAIVE_EST(dt), error = function(e) NULL)
  if (is.null(ans)) return(.na_result("MAIVE model failed."))

  co <- summary(ans$model)$coefficients
  ci <- tryCatch(stats::confint(ans$model), error = function(e) matrix(NA_real_, nrow = 2, ncol = 2))
  est <- .make_estimates(
    term = c("b0", "se"),
    estimate = as.numeric(co[, 1]),
    std_error = as.numeric(co[, 2]),
    statistic = as.numeric(co[, 3]),
    p_value = as.numeric(co[, 4]),
    conf_low = as.numeric(ci[, 1]),
    conf_high = as.numeric(ci[, 2])
  )
  list(estimates = est, pub_bias_p = est$p_value[est$term == "se"], status = "ok", note = "")
}

.fit_fativ <- function(dt) {
  if (is.null(dt$obs)) {
    return(.na_result("EGGER5 / FATIV requires `obs`.", status = "requires_obs"))
  }
  ans <- tryCatch(.MAIVE2_EST(dt), error = function(e) NULL)
  if (is.null(ans)) return(.na_result("FATIV model failed."))

  co <- summary(ans$model)$coefficients
  ci <- tryCatch(stats::confint(ans$model), error = function(e) matrix(NA_real_, nrow = 2, ncol = 2))
  est <- .make_estimates(
    term = c("b0", "se"),
    estimate = as.numeric(co[, 1]),
    std_error = as.numeric(co[, 2]),
    statistic = as.numeric(co[, 3]),
    p_value = as.numeric(co[, 4]),
    conf_low = as.numeric(ci[, 1]),
    conf_high = as.numeric(ci[, 2])
  )
  list(estimates = est, pub_bias_p = est$p_value[est$term == "se"], status = "ok", note = "")
}

# Selection models ---------------------------------------------------------

.fit_psm3 <- function(dt) {
  psm3 <- tryCatch({
    myRE <- metafor::rma(yi = dt$eff, vi = dt$se^2, method = "REML")
    metafor::selmodel(myRE, type = "stepfun", steps = c(0.025, 1))
  }, error = function(e) NULL)
  if (is.null(psm3)) return(.na_result("3PSM failed or did not converge."))

  cv.025 <- stats::qnorm(0.975) * dt$se
  est <- c(psm3$tau2[1], psm3$beta[1], psm3$delta[2])
  se <- c(psm3$se.tau2, psm3$se, psm3$se.delta[2])
  pub.sel.3psm <- mean(
    stats::pnorm(cv.025, est[2], sqrt(est[1] + dt$se^2), lower.tail = FALSE) +
      est[3] * stats::pnorm(cv.025, est[2], sqrt(est[1] + dt$se^2), lower.tail = TRUE)
  )
  out <- .make_estimates(
    term = c("tau2", "b0", "pr.nonsig", "pubSelRate"),
    estimate = c(est, pub.sel.3psm),
    std_error = c(se, NA_real_),
    statistic = c(est / se, NA_real_),
    p_value = c(stats::pnorm(est / se, lower.tail = FALSE) * 2, NA_real_),
    conf_low = c(est + stats::qnorm(0.025) * se, NA_real_),
    conf_high = c(est + stats::qnorm(0.975) * se, NA_real_)
  )
  list(estimates = out, pub_bias_p = psm3$LRTp, status = "ok", note = "")
}

.fit_psm4 <- function(dt) {
  psm4 <- tryCatch({
    myRE <- metafor::rma(yi = dt$eff, vi = dt$se^2, method = "REML")
    metafor::selmodel(myRE, type = "stepfun", steps = c(0.025, 0.5, 1))
  }, error = function(e) NULL)
  if (is.null(psm4)) return(.na_result("4PSM failed or did not converge."))

  cv.025 <- stats::qnorm(0.975) * dt$se
  est <- c(psm4$tau2[1], psm4$beta[1], psm4$delta[2], psm4$delta[3])
  se <- c(psm4$se.tau2, psm4$se, psm4$se.delta[2], psm4$se.delta[3])

  sig.plus <- stats::pnorm(cv.025, est[2], sqrt(est[1] + dt$se^2), lower.tail = FALSE)
  nonsig.plus <- stats::pnorm(0, est[2], sqrt(est[1] + dt$se^2), lower.tail = FALSE) - sig.plus
  minus <- stats::pnorm(0, est[2], sqrt(est[1] + dt$se^2), lower.tail = TRUE)
  pub.sel.4psm <- mean(sig.plus + est[3] * nonsig.plus + est[4] * minus)

  out <- .make_estimates(
    term = c("tau2", "b0", "pr.nonsig", "pr.opposite", "pubSelRate"),
    estimate = c(est, pub.sel.4psm),
    std_error = c(se, NA_real_),
    statistic = c(est / se, NA_real_),
    p_value = c(stats::pnorm(est / se, lower.tail = FALSE) * 2, NA_real_),
    conf_low = c(est + stats::qnorm(0.025) * se, NA_real_),
    conf_high = c(est + stats::qnorm(0.975) * se, NA_real_)
  )
  list(estimates = out, pub_bias_p = psm4$LRTp, status = "ok", note = "")
}

# AK models ----------------------------------------------------------------

.AK1logLik <- function(para, mydata, z = c("est", "log")) {
  z <- match.arg(z)
  n <- nrow(mydata)
  tauhat <- para[1]
  betap <- para[2]
  Coeff <- as.matrix(para[3:length(para)])
  if (ncol(mydata) > 3) {
    err <- mydata[, 2] - as.matrix(mydata[, c(4:ncol(mydata))]) %*% Coeff
    esthat <- as.matrix(mydata[, c(4:ncol(mydata))]) %*% Coeff
  } else {
    err <- mydata[, 2] - Coeff
    esthat <- Coeff
  }

  se <- mydata[, 3]
  t <- mydata[, 2] / mydata[, 3]
  cutoffs <- c(-1.96, 1.96)
  phat <- (abs(t) < 1.96) * betap + (abs(t) >= 1.96) * 1
  meanbeta <- rep(0, n)
  for (i in seq_len(n)) {
    prob_mid <- (stats::pnorm((cutoffs[2] * se[i] - esthat[i]) / sqrt((tauhat) + se[i]^2)) -
                   stats::pnorm((cutoffs[1] * se[i] - esthat[i]) / sqrt((tauhat) + se[i]^2)))
    prob_ex <- 1 - prob_mid
    meanbeta[i] <- betap * prob_mid + 1 * prob_ex
  }
  fX <- stats::dnorm(err, 0, sqrt(se^2 + tauhat))
  L <- (phat / meanbeta) * fX
  logL <- log(L)
  LLH <- -sum(log(L))
  if (tauhat < 0 || betap < 0) {
    LLH <- 1e10
  }
  if (z == "est") LLH else logL
}

.AK2logLik <- function(para, mydata, z = c("est", "log")) {
  z <- match.arg(z)
  n <- nrow(mydata)
  tauhat <- para[1]
  beta1 <- para[2]
  beta2 <- para[3]
  beta3 <- para[4]
  Coeff <- as.matrix(para[5:length(para)])
  if (ncol(mydata) > 3) {
    err <- mydata[, 2] - as.matrix(mydata[, c(4:ncol(mydata))]) %*% Coeff
    esthat <- as.matrix(mydata[, c(4:ncol(mydata))]) %*% Coeff
  } else {
    err <- mydata[, 2] - Coeff
    esthat <- Coeff
  }

  se <- mydata[, 3]
  t <- mydata[, 2] / mydata[, 3]
  cutoffs <- c(-1.96, 0, 1.96)
  phat <- (t <= -1.96) * beta1 + (t > -1.96 & t <= 0) * beta2 + (0 < t & t < 1.96) * beta3 + (t >= 1.96) * 1
  meanbeta <- rep(0, n)
  for (i in seq_len(n)) {
    prob_vlow <- stats::pnorm((cutoffs[1] * se[i] - esthat[i]) / sqrt((tauhat) + se[i]^2))
    prob_low <- (stats::pnorm((cutoffs[2] * se[i] - esthat[i]) / sqrt((tauhat) + se[i]^2)) -
                   stats::pnorm((cutoffs[1] * se[i] - esthat[i]) / sqrt((tauhat) + se[i]^2)))
    prob_upper <- (stats::pnorm((cutoffs[3] * se[i] - esthat[i]) / sqrt((tauhat) + se[i]^2)) -
                     stats::pnorm((cutoffs[2] * se[i] - esthat[i]) / sqrt((tauhat) + se[i]^2)))
    prob_vupper <- 1 - prob_vlow - prob_low - prob_upper
    meanbeta[i] <- beta1 * prob_vlow + beta2 * prob_low + beta3 * prob_upper + 1 * prob_vupper
  }

  fX <- stats::dnorm(err, 0, sqrt(se^2 + tauhat))
  L <- (phat / meanbeta) * fX
  logL <- log(L)
  LLH <- -sum(log(L))
  if (tauhat < 0 || beta1 < 0 || beta2 < 0 || beta3 < 0) {
    LLH <- 1e10
  }
  if (z == "est") LLH else logL
}

.fit_ak1 <- function(dt) {
  AKdata <- as.data.frame(cbind(id = as.numeric(dt$StdID), effect = dt$eff, se = dt$se, constant = 1))
  InitialValue <- as.numeric(mean(AKdata$effect))
  t2 <- mean(AKdata$se^2) / 4
  fn <- function(par) .AK1logLik(par, AKdata, "est")
  Result <- tryCatch(
    nlminb(fn,
           start = c(t2, 0, InitialValue),
           lower = c(0, 0.01, rep(-Inf, length(InitialValue))),
           upper = c(Inf, Inf, rep(Inf, length(InitialValue))),
           hessian = TRUE,
           control = list(iter.max = 1000, abs.tol = 10^(-20), eval.max = 1000)),
    error = function(e) NULL
  )
  if (is.null(Result) || Result$convergence == 1) {
    return(.na_result("AK1 failed or did not converge."))
  }

  se_vec <- tryCatch(
    diag(sqrt(MASS::ginv(numDeriv::hessian(fn, Result$par), tol = 10^(-30)))),
    error = function(e) rep(NA_real_, length(Result$par))
  )
  estcoef <- as.numeric(Result$par)
  if (any(!is.finite(se_vec))) {
    return(.na_result("AK1 Hessian inversion failed."))
  }
  crit <- stats::qt(0.975, df = (nrow(AKdata) - length(estcoef)))
  b0_est <- estcoef[3]
  b0_se <- se_vec[3]
  rho1_est <- estcoef[2]
  rho1_se <- se_vec[2]
  tau2_est <- estcoef[1]
  tau2_se <- se_vec[1]
  AK1LLK <- .AK1logLik(Result$par, AKdata, "est")

  re_mod <- tryCatch(metafor::rma(dt$eff, dt$se^2, method = "ML"), error = function(e) NULL)
  RE_LLK <- if (is.null(re_mod)) NA_real_ else as.numeric(logLik(re_mod, REML = FALSE))
  cv.025 <- stats::qnorm(0.975) * dt$se
  if (!is.na(RE_LLK) && !is.na(AK1LLK)) {
    LRT <- -2 * (AK1LLK - RE_LLK)
    pLRT <- stats::pchisq(LRT, 1, lower.tail = FALSE)
    sig.plus <- stats::pnorm(cv.025, b0_est, sqrt(tau2_est + dt$se^2), lower.tail = FALSE)
    nonsig.plus <- stats::pnorm(0, b0_est, sqrt(tau2_est + dt$se^2), lower.tail = FALSE) - sig.plus
    sig.min <- stats::pnorm(-cv.025, b0_est, sqrt(tau2_est + dt$se^2), lower.tail = TRUE)
    nonsig.min <- stats::pnorm(0, b0_est, sqrt(tau2_est + dt$se^2), lower.tail = TRUE) - sig.min
    pub.sel <- mean(sig.plus + rho1_est * (nonsig.plus + nonsig.min) + sig.min)
  } else {
    pLRT <- NA_real_
    pub.sel <- NA_real_
  }

  est <- rbind(
    .make_estimates("b0", b0_est, b0_se, b0_est / b0_se,
                    stats::dt(b0_est / b0_se, df = nrow(AKdata) - length(estcoef)),
                    b0_est - crit * b0_se, b0_est + crit * b0_se),
    .make_estimates("rho1", rho1_est, rho1_se),
    .make_estimates("tau2", tau2_est, tau2_se),
    .make_estimates("pubSelRate", pub.sel)
  )
  list(estimates = est, pub_bias_p = pLRT, status = "ok", note = "")
}

.fit_ak2 <- function(dt) {
  AKdata <- as.data.frame(cbind(id = as.numeric(dt$StdID), effect = dt$eff, se = dt$se, constant = 1))
  InitialValue <- as.numeric(mean(AKdata$effect))
  t2 <- mean(AKdata$se^2) / 4
  fn <- function(par) .AK2logLik(par, AKdata, "est")
  Result <- tryCatch(
    nlminb(fn,
           start = c(t2, 0, 0, 0, InitialValue),
           lower = c(0, 0.01, 0.01, 0.01, rep(-Inf, length(InitialValue))),
           upper = c(Inf, Inf, Inf, Inf, rep(Inf, length(InitialValue))),
           hessian = TRUE,
           control = list(iter.max = 1000, abs.tol = 10^(-20), eval.max = 1000)),
    error = function(e) NULL
  )
  if (is.null(Result) || Result$convergence == 1) {
    return(.na_result("AK2 failed or did not converge."))
  }

  se_vec <- tryCatch(
    diag(sqrt(MASS::ginv(numDeriv::hessian(fn, Result$par), tol = 10^(-30)))),
    error = function(e) rep(NA_real_, length(Result$par))
  )
  estcoef <- as.numeric(Result$par)
  if (any(!is.finite(se_vec))) {
    return(.na_result("AK2 Hessian inversion failed."))
  }

  crit <- stats::qt(0.975, df = (nrow(AKdata) - length(estcoef)))
  tau2_est <- estcoef[1]; tau2_se <- se_vec[1]
  rho1_est <- estcoef[2]; rho1_se <- se_vec[2]
  rho2_est <- estcoef[3]; rho2_se <- se_vec[3]
  rho3_est <- estcoef[4]; rho3_se <- se_vec[4]
  b0_est <- estcoef[5]; b0_se <- se_vec[5]
  AK2LLK <- .AK2logLik(Result$par, AKdata, "est")

  re_mod <- tryCatch(metafor::rma(dt$eff, dt$se^2, method = "ML"), error = function(e) NULL)
  RE_LLK <- if (is.null(re_mod)) NA_real_ else as.numeric(logLik(re_mod, REML = FALSE))
  cv.025 <- stats::qnorm(0.975) * dt$se
  if (!is.na(RE_LLK) && !is.na(AK2LLK)) {
    pLRT <- stats::pchisq(-2 * (AK2LLK - RE_LLK), 2, lower.tail = FALSE)
    sig.plus <- stats::pnorm(cv.025, b0_est, sqrt(tau2_est + dt$se^2), lower.tail = FALSE)
    nonsig.plus <- stats::pnorm(0, b0_est, sqrt(tau2_est + dt$se^2), lower.tail = FALSE) - sig.plus
    sig.min <- stats::pnorm(-cv.025, b0_est, sqrt(tau2_est + dt$se^2), lower.tail = TRUE)
    nonsig.min <- stats::pnorm(0, b0_est, sqrt(tau2_est + dt$se^2), lower.tail = TRUE) - sig.min
    pub.sel <- mean(sig.plus + rho3_est * nonsig.plus + rho2_est * nonsig.min + rho1_est * sig.min)
  } else {
    pLRT <- NA_real_
    pub.sel <- NA_real_
  }

  est <- rbind(
    .make_estimates("b0", b0_est, b0_se, b0_est / b0_se,
                    stats::dt(b0_est / b0_se, df = nrow(AKdata) - length(estcoef)),
                    b0_est - crit * b0_se, b0_est + crit * b0_se),
    .make_estimates("rho1", rho1_est, rho1_se),
    .make_estimates("rho2", rho2_est, rho2_se),
    .make_estimates("rho3", rho3_est, rho3_se),
    .make_estimates("tau2", tau2_est, tau2_se),
    .make_estimates("pubSelRate", pub.sel)
  )
  list(estimates = est, pub_bias_p = pLRT, status = "ok", note = "")
}

# EK, TES, PSST/TESS, BEGG, p-uniform -------------------------------------

.fit_ek <- function(dt, tb = 1.96, ts = 1.96) {
  d <- dt$eff
  v <- dt$se^2
  se <- sqrt(v)
  m <- length(v)
  out <- tryCatch({
    FP1 <- summary(stats::lm(d ~ se, weights = 1 / v))
    if (FP1$coefficient[1, 4] > 0.05) {
      alpha1 <- as.numeric(FP1$coefficient[1, 1])
      Q <- sum((as.numeric(FP1$residuals) / se)^2)
    } else {
      FP2 <- summary(stats::lm(d ~ v, weights = 1 / v))
      alpha1 <- as.numeric(FP2$coefficient[1, 1])
      Q <- sum((as.numeric(FP2$residuals) / se)^2)
    }
    wi <- 1 / v
    sig_eta2 <- max(0, m * ((Q / (m - 2)) - 1) / sum(wi))
    a <- (((alpha1^2) - tb * tb * sig_eta2) / ((tb + ts) * alpha1)) * (alpha1 > (tb * sqrt(sig_eta2)))
    g <- (se - a) * (se >= a)
    reg <- stats::lm(d ~ g, weights = 1 / v)
    EKreg <- summary(reg)$coefficient
    EKConf <- stats::confint(reg)
    list(EKreg = EKreg, EKConf = EKConf)
  }, error = function(e) NULL)

  if (is.null(out) || any(is.na(out$EKConf))) return(.na_result("EK failed."))

  est <- .make_estimates(
    term = c("b0", "g"),
    estimate = as.numeric(out$EKreg[, 1]),
    std_error = as.numeric(out$EKreg[, 2]),
    statistic = as.numeric(out$EKreg[, 3]),
    p_value = as.numeric(out$EKreg[, 4]),
    conf_low = as.numeric(out$EKConf[, 1]),
    conf_high = as.numeric(out$EKConf[, 2])
  )
  list(estimates = est, pub_bias_p = est$p_value[est$term == "g"], status = "ok", note = "")
}

.fit_tes <- function(dt, alpha = 0.05, side = "right") {
  out <- tryCatch({
    est.fe <- metafor::rma(yi = dt$eff, sei = dt$se, method = "FE")$b[1]
    if (side == "right") {
      pow <- stats::pnorm(stats::qnorm(alpha, lower.tail = FALSE, sd = dt$se), mean = est.fe,
                          sd = dt$se, lower.tail = FALSE)
      O <- sum(stats::pnorm(dt$eff / dt$se, lower.tail = FALSE) < alpha)
    } else {
      pow <- stats::pnorm(stats::qnorm(alpha, sd = dt$se), mean = est.fe, sd = dt$se)
      O <- sum(stats::pnorm(dt$eff / dt$se) < alpha)
    }
    E <- sum(pow)
    n <- length(dt$eff)
    A <- (O - E)^2 / E + (O - E)^2 / (n - E)
    pval.tes <- stats::pchisq(A, 1, lower.tail = FALSE)
    pval.tes <- ifelse(pval.tes < 0.5, pval.tes * 2, (1 - pval.tes) * 2)
    list(A = A, pval.tes = pval.tes, O = O, E = E, n = n)
  }, error = function(e) NULL)

  if (is.null(out)) return(.na_result("TES failed."))
  list(estimates = NULL, pub_bias_p = out$pval.tes, status = "ok", note = "")
}

.fit_psst_tess <- function(dt) {
  out <- tryCatch({
    k <- length(dt$eff)
    t <- dt$eff / dt$se
    Precision <- 1 / dt$se
    MA <- metafor::rma(yi = dt$eff, sei = dt$se, method = "DL")
    HetVar <- MA$tau2
    reg <- stats::lm(t ~ 0 + Precision)
    UWLS <- as.numeric(reg$coefficients)
    zz <- ((1.96 * dt$se - UWLS) / (dt$se^2 + HetVar)^0.5)
    Esigtot <- sum(1 - stats::pnorm(zz))
    SS <- (t > 1.96) * 1
    SStot <- sum(SS)
    Pss <- SStot / k
    ESS <- (SStot - Esigtot) / k
    Pe <- Esigtot / k
    PSST <- (Pss - Pe) / (Pe * (1 - Pe) / k)^0.5
    TESS <- (ESS - 0.05) / (0.0475 / k)^0.5
    list(
      psst_p = (1 - stats::pnorm(PSST)),
      tess_p = (1 - stats::pnorm(TESS)),
      psst_z = PSST,
      tess_z = TESS
    )
  }, error = function(e) NULL)

  if (is.null(out)) {
    return(list(psst_p = NA_real_, tess_p = NA_real_, status = "failed", note = "PSST/TESS failed."))
  }
  list(psst_p = out$psst_p, tess_p = out$tess_p, status = "ok", note = "")
}

.fit_begg <- function(dt) {
  out <- tryCatch(
    metafor::ranktest(dt$eff, sei = dt$se, exact = FALSE)$pval,
    error = function(e) NA_real_
  )
  list(pub_bias_p = as.numeric(out), status = ifelse(is.na(out), "failed", "ok"), note = ifelse(is.na(out), "Begg test failed.", ""))
}

.fit_puniform <- function(dt) {
  PuniTest <- tryCatch(
    puniform::puniform(yi = dt$eff, vi = dt$se^2, side = "right"),
    error = function(e) NULL
  )
  if (is.null(PuniTest)) return(.na_result("p-uniform failed."))
  est <- .make_estimates(
    term = "b0",
    estimate = as.numeric(PuniTest$est),
    std_error = NA_real_,
    statistic = NA_real_,
    p_value = as.numeric(PuniTest$pval.0),
    conf_low = as.numeric(PuniTest$ci.lb),
    conf_high = as.numeric(PuniTest$ci.ub)
  )
  list(estimates = est, pub_bias_p = as.numeric(PuniTest$pval.pb), status = "ok", note = "")
}

# Skewness ----------------------------------------------------------------

.fit_skew <- function(dt, reMA_model) {
  if (is.null(reMA_model)) {
    return(list(skew_p = NA_real_, combined_p = NA_real_, status = "failed", note = "Baseline RE model unavailable."))
  }
  out <- tryCatch({
    tau2 <- reMA_model$tau2
    if (!is.finite(tau2)) stop("tau2 is not finite")
    std.y <- dt$eff / sqrt(dt$se^2 + tau2)
    std.x <- 1 / sqrt(dt$se^2 + tau2)
    skew.reg <- stats::lm(std.y ~ std.x)
    std.dev <- as.numeric(summary(skew.reg)$residuals)
    cm2 <- stats::var(std.dev)
    cm3 <- mean((std.dev - mean(std.dev))^3)
    cm4 <- mean((std.dev - mean(std.dev))^4)
    cm5 <- mean((std.dev - mean(std.dev))^5)
    cm6 <- mean((std.dev - mean(std.dev))^6)
    skewness <- cm3 / (cm2^(1.5))
    var0.skew <- 6
    skewness.pval <- 2 * (1 - stats::pnorm(sqrt(length(dt$eff) / var0.skew) * abs(skewness)))
    reg.coef <- summary(skew.reg)$coef
    reg.pval <- reg.coef["(Intercept)", "Pr(>|t|)"]
    combined.pval <- 1 - (1 - min(c(reg.pval, skewness.pval)))^2
    list(skewness.pval = skewness.pval, combined.pval = combined.pval)
  }, error = function(e) NULL)

  if (is.null(out)) {
    return(list(skew_p = NA_real_, combined_p = NA_real_, status = "failed", note = "Skewness tests failed."))
  }
  list(skew_p = out$skewness.pval, combined_p = out$combined.pval, status = "ok", note = "")
}

# Calipers ----------------------------------------------------------------

.fit_calipers <- function(dt, alpha = 0.05) {
  zcrit <- stats::qnorm(1 - alpha / 2)
  abtstat <- abs(dt$eff / dt$se)
  windows <- c(0.05, 0.10, 0.15, 0.20)
  code_methods <- c("Caliper05", "Caliper10", "Caliper15", "Caliper20")
  notes <- character(length(windows))
  pvals <- numeric(length(windows))
  status <- character(length(windows))

  for (i in seq_along(windows)) {
    w <- windows[i]
    over <- sum(zcrit < abtstat & abtstat < zcrit * (1 + w))
    under <- sum(zcrit * (1 - w) < abtstat & abtstat < zcrit)
    total <- over + under
    if (total < 5) {
      pvals[i] <- NA_real_
      status[i] <- "not_applicable"
      notes[i] <- "Fewer than 5 studies in the combined caliper interval."
    } else {
      pvals[i] <- stats::binom.test(c(over, under), alternative = "greater")$p.value
      status[i] <- "ok"
      notes[i] <- paste0("Counts: over=", over, ", under=", under)
    }
  }

  list(summary = data.frame(
    code_method = code_methods,
    publication_bias_p = pvals,
    status = status,
    note = notes,
    stringsAsFactors = FALSE
  ))
}
