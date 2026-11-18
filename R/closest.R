#' Find individual simulation conditions closest to a meta-analysis
#'
#' Filters the full simulation conditions, matches K first, then brackets the
#' target tau using each condition's mean estimated tau. Censoring and QRP
#' levels remain separate throughout selection. Performance is summarised only
#' after the neighbours have been selected.
#'
#' @param k A single finite whole number of studies, at least 3.
#' @param tau Estimated between-study standard deviation, finite and non-negative.
#'   Supply `sqrt(tau2)` if your estimate is a variance.
#' @param design `"extended"` (864 conditions, K = 5 to 800) or `"original"`
#'   (432 conditions, K = 10 to 100, based on Carter et al., 2019).
#' @param qrp `"all"`, `"without"` (none), or `"with"` (medium/high).
#' @param universe Legacy additional filter: `"all"` or `"pb_only"` (QRP or
#'   Censoring medium/high). The standard app uses the direct QRP/Censoring filters.
#' @param pooling Retained to give existing scripts a clear migration error.
#'   Only `"none"` is supported. Pooling before matching is no longer performed.
#' @param censoring `"all"`, `"without"` (none), or `"with"` (medium/high),
#'   independently of QRP. Selecting All keeps all three Censoring levels separate.
#'
#' @details Each condition has 1,000 simulation replications. `mean_tau` is the
#'   mean of valid estimated standard deviations within that condition;
#'   `mean_tau2` is the mean estimated variance, not squared mean tau.
#'   No Censoring levels, QRP levels or conditions are pooled before matching.
#'   Exact K is used when available; otherwise the nearest lower and upper K
#'   are searched separately. At each K the immediately lower/upper mean tau
#'   is returned. Exact tau matches take priority; ties within 1e-12 on the tau
#'   scale are retained. Only available boundaries are returned outside the range.
#'
#'   A condition is biased when QRP or Censoring is medium/high. Individual
#'   biased conditions provide TPR; individual unbiased conditions provide FPR.
#'   An individual condition's Balanced Accuracy is unavailable.
#'
#'   After selection, `summary$mean_tpr` averages the selected biased conditions'
#'   rejection rates, and `mean_fpr` averages selected unbiased conditions'
#'   rejection rates. `balanced_accuracy = (mean_tpr + 1 - mean_fpr) / 2`, giving
#'   the two classes equal weight. Each condition receives equal weight within
#'   its class. Both classes are required; all selected conditions within each
#'   class must have valid rates. Missing results are never replaced with zero
#'   and no extra control conditions are introduced. This is a descriptive
#'   comparison of the selected neighbours, not an interpolated prediction.
#'   Mean rejection rate and applicability give equal weight to all matches.
#'   Low applicability is reported separately.
#'
#' @return A `pbias_closest_result` with `matches` (full individual conditions,
#'   including `bias`, an explicit QRP/Censoring label such as
#'   `QRP = None & Censoring = High`),
#'   `performance` (21 numeric method results per match), `summary` (method
#'   comparison after selection, sorted by Balanced Accuracy), `candidates`
#'   (eligible individual conditions), `notes`, and `meta` (inputs, searched K
#'   counts, weighting and simulation provenance). Rates are proportions.
#' @export
#' @examples
#' near <- pbias_closest(k = 65, tau = 0.2)
#' near
#' near$summary
#' near$matches
#' pbias_closest(60, sqrt(0.04), qrp = "without", censoring = "all")
pbias_closest <- function(k, tau, design = c("extended", "original"),
                          qrp = c("all", "without", "with"),
                          universe = c("all", "pb_only"),
                          pooling = c("none", "publication_bias"),
                          censoring = c("all", "without", "with")) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 3 || k != floor(k)) {
    stop("k must be one finite whole number of studies, at least 3.", call. = FALSE)
  }
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0) {
    stop("tau must be one finite, non-negative estimated standard deviation; use sqrt(tau2) for a variance.", call. = FALSE)
  }
  design <- match.arg(design); qrp <- match.arg(qrp)
  universe <- match.arg(universe); pooling <- match.arg(pooling)
  censoring <- match.arg(censoring)
  if (pooling != "none") stop("Censoring pooling before matching is no longer supported. Omit pooling or use pooling = 'none'.", call. = FALSE)
  data <- readRDS(.pbias_extdata_file("closest_conditions.rds"))
  candidates <- data$conditions
  bias_levels <- c(none = "None", med = "Medium", high = "High")
  candidates$bias <- paste0("QRP = ", unname(bias_levels[candidates$QRP]),
                            " & Censoring = ", unname(bias_levels[candidates$Censor]))
  candidates$profile_id <- paste0("condition-", candidates$Environment)
  candidates$n_conditions <- 1L
  candidates$n_unbiased_conditions <- as.integer(!candidates$biased)
  candidates$n_biased_conditions <- as.integer(candidates$biased)
  if (universe == "pb_only") candidates <- candidates[candidates$biased, , drop = FALSE]
  if (censoring == "without") candidates <- candidates[candidates$Censor == "none", , drop = FALSE]
  if (censoring == "with") candidates <- candidates[candidates$Censor != "none", , drop = FALSE]
  if (design == "original") candidates <- candidates[candidates$K >= 10 & candidates$K <= 100, , drop = FALSE]
  if (qrp == "without") candidates <- candidates[candidates$QRP == "none", , drop = FALSE]
  if (qrp == "with") candidates <- candidates[candidates$QRP != "none", , drop = FALSE]
  if (!nrow(candidates)) stop("No simulation conditions match these filters. Use universe = 'all' to include unbiased conditions.", call. = FALSE)
  candidates <- candidates[order(candidates$K, candidates$mean_tau, candidates$profile_id), , drop = FALSE]
  rownames(candidates) <- NULL
  matches <- .pbias_bracket(candidates, k, tau)
  matches$match_id <- paste0("M", seq_len(nrow(matches)))
  matches$match_label <- paste0(matches$match_id, ": K = ", matches$K,
                              ", mean tau = ", formatC(matches$mean_tau, format = "f", digits = 6))
  matches$role <- ifelse(matches$biased, "Biased", "Unbiased")
  perf <- data$condition_performance
  perf <- perf[perf$Environment %in% matches$Environment, , drop = FALSE]
  index <- match(perf$Environment, matches$Environment)
  perf$profile_id <- matches$profile_id[index]
  perf$match_id <- matches$match_id[index]
  perf$n_total <- matches$n_total[index]
  perf$rejection_rate <- perf$rate
  perf$applicability <- perf$n_valid / perf$n_total
  perf$tpr <- ifelse(matches$biased[index], perf$rate, NA_real_)
  perf$fpr <- ifelse(!matches$biased[index], perf$rate, NA_real_)
  perf$balanced_accuracy <- NA_real_
  perf$status <- ifelse(perf$n_valid == 0, "Unavailable",
                       ifelse(perf$applicability < 0.9, "Below 90% applicability", "Meets 90% threshold"))
  methods <- c("BEGG", paste0("EGGER", 1:5), "EK", "SKEW1", "SKEW2", "PSM3", "PSM4",
               "AK1", "AK2", "PUNIF", "TES", "PSST", "TESS", "CALI05", "CALI10", "CALI15", "CALI20")
  perf <- perf[order(match(perf$match_id, matches$match_id), match(perf$method, methods)),
               c("match_id", "profile_id", "method", "rejection_rate", "applicability", "n_valid", "n_total",
                 "tpr", "fpr", "balanced_accuracy", "status"), drop = FALSE]
  rownames(perf) <- NULL
  summary <- .pbias_closest_summary(perf, matches$match_id, matches$biased)
  ks <- sort(unique(candidates$K)); searched_k <- sort(unique(matches$K))
  searched_counts <- table(factor(candidates$K[candidates$K %in% searched_k], levels = searched_k))
  notes <- character()
  if (k < min(ks) || k > max(ks)) notes <- c(notes, paste0("K is outside the simulated range (", min(ks), "-", max(ks), "). Only the nearest available K is shown."))
  for (value_k in searched_k) {
    values <- candidates$mean_tau[candidates$K == value_k]
    if (tau < min(values) || tau > max(values)) notes <- c(notes, paste0("At K = ", value_k, ", tau is outside the available range; only the nearest boundary condition is shown."))
  }
  structure(list(matches = matches, performance = perf, summary = summary, candidates = candidates, notes = notes,
    meta = c(list(k = k, tau = tau, design = design, qrp = qrp, universe = universe, pooling = "none",
                  censoring = censoring, package_version = "0.4.1", n_candidates = nrow(candidates),
                  searched_k = searched_k, searched_counts = stats::setNames(as.integer(searched_counts), searched_k),
                  n_searched = sum(searched_counts),
                  summary_weighting = "Mean TPR across selected biased conditions; mean FPR across selected unbiased conditions; equal weight per class; complete rates required",
                  matching_unit = "Individual full simulation condition; no pooling before selection"), data$metadata)),
    class = "pbias_closest_result")
}

.pbias_bracket <- function(candidates, k, tau) {
  ks <- sort(unique(candidates$K))
  selected_k <- if (k %in% ks) k else c(if (any(ks < k)) max(ks[ks < k]),
                                        if (any(ks > k)) min(ks[ks > k]))
  tolerance <- 1e-12 * max(1, abs(tau))
  result <- lapply(selected_k, function(value_k) {
    p <- candidates[candidates$K == value_k, , drop = FALSE]
    exact <- abs(p$mean_tau - tau) <= tolerance
    if (any(exact)) {
      z <- p[exact, , drop = FALSE]
      z$heterogeneity_match <- "Exact tau"
    } else {
      lower <- p$mean_tau[p$mean_tau < tau]
      upper <- p$mean_tau[p$mean_tau > tau]
      lo <- if (length(lower)) p[abs(p$mean_tau - max(lower)) <= tolerance, , drop = FALSE] else NULL
      hi <- if (length(upper)) p[abs(p$mean_tau - min(upper)) <= tolerance, , drop = FALSE] else NULL
      if (!is.null(lo)) lo$heterogeneity_match <- "Tau below"
      if (!is.null(hi)) hi$heterogeneity_match <- "Tau above"
      z <- rbind(lo, hi)
    }
    z$k_match <- if (value_k == k) "Exact K" else if (value_k < k) "K below" else "K above"
    z
  })
  out <- do.call(rbind, result)
  rownames(out) <- NULL
  out
}

#' Print closest simulation conditions
#' @param x An object returned by [pbias_closest()].
#' @param ... Further arguments (unused).
#' @return `x`, invisibly.
#' @export
print.pbias_closest_result <- function(x, ...) {
  cat("Closest individual simulation conditions\n")
  cat(sprintf("K = %s; estimated tau = %.6g; tau-squared = %.6g\n", x$meta$k, x$meta$tau, x$meta$tau^2))
  cat(sprintf("%s design | QRP: %s | Censoring: %s\n", x$meta$design, x$meta$qrp, x$meta$censoring))
  cat(sprintf("%s eligible conditions; %s searched at K = %s.\n", x$meta$n_candidates, x$meta$n_searched, paste(x$meta$searched_k, collapse = ", ")))
  m <- x$matches
  print(data.frame(Match = m$match_id, K = m$K, Mean_tau = round(m$mean_tau, 6), Mean_tau2 = round(m$mean_tau2, 6),
                   Bias = m$bias), row.names = FALSE)
  s <- x$summary
  pct <- function(z) ifelse(is.finite(z), sprintf("%.1f%%", 100*z), "n/a")
  cat("\nDescriptive performance across the selected conditions:\n")
  print(data.frame(Method = s$method, BA = pct(s$balanced_accuracy), TPR = pct(s$mean_tpr), FPR = pct(s$mean_fpr),
                   TPR_matches = paste0(s$tpr_available, "/", s$n_biased_matches),
                   FPR_matches = paste0(s$fpr_available, "/", s$n_unbiased_matches)), row.names = FALSE)
  cat("BA = (mean TPR + 1 - mean FPR) / 2; both classes and complete selected rates required.\n")
  cat("Per-condition BA is unavailable. Inspect $performance for rates/applicability and $matches for settings.\n")
  if (length(x$notes)) cat(paste(x$notes, collapse = "\n"), "\n")
  invisible(x)
}

.pbias_closest_summary <- function(performance, match_ids, biased) {
  methods <- unique(performance$method); n <- length(match_ids)
  complete_mean <- function(z) if (length(z) && all(is.finite(z))) mean(z) else NA_real_
  out <- do.call(rbind, lapply(methods, function(method) {
    p <- performance[performance$method == method, , drop = FALSE]
    p <- p[match(match_ids, p$match_id), , drop = FALSE]
    tpr <- complete_mean(p$rejection_rate[biased]); fpr <- complete_mean(p$rejection_rate[!biased])
    ba <- (tpr + 1 - fpr)/2
    data.frame(method = method, balanced_accuracy = ba, mean_tpr = tpr, mean_fpr = fpr,
      mean_rejection_rate = complete_mean(p$rejection_rate), mean_applicability = complete_mean(p$applicability),
      min_applicability = if (all(is.finite(p$applicability))) min(p$applicability) else NA_real_,
      n_matches = n, n_biased_matches = sum(biased), n_unbiased_matches = sum(!biased),
      tpr_available = sum(is.finite(p$rejection_rate[biased])), fpr_available = sum(is.finite(p$rejection_rate[!biased])),
      rejection_available = sum(is.finite(p$rejection_rate)),
      low_applicability_matches = sum(!is.finite(p$applicability) | p$applicability < .9),
      comparison_status = if (!any(biased) || all(biased)) "Needs biased and unbiased matches" else
        if (!is.finite(ba)) "Missing result in a selected condition" else "Complete selected-condition comparison")
  }))
  out <- out[order(-out$balanced_accuracy, match(out$method, methods), na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}
