#' Find simulation conditions closest to a meta-analysis
#'
#' Matches the number of studies first, then brackets the supplied estimated
#' heterogeneity with mean estimated tau from the simulations. Results for
#' neighbouring conditions are kept separate; no distance score or interpolation
#' is used. This is a lookup of simulation performance, not a test of new data.
#'
#' @param k A single whole number of studies, at least 3.
#' @param tau A single finite, non-negative estimated between-study standard
#'   deviation. If you have estimated tau-squared, supply `sqrt(tau2)`.
#' @param design `"extended"` (864 conditions, K = 5 to 800) or `"original"`
#'   (432 conditions, K = 10 to 100, based on Carter et al., 2019).
#' @param qrp `"without"` for no p-hacking, `"with"` for medium/high p-hacking,
#'   or `"all"` for all three QRP levels. Levels remain separate.
#' @param universe `"all"` (default) includes all conditions. `"pb_only"`
#'   includes any biased condition: medium/high censoring OR medium/high QRP,
#'   including QRP-only conditions without censoring.
#' @param pooling `"publication_bias"` (default) holds K, true effect, population
#'   tau and QRP fixed and averages eligible censoring levels equally. All uses
#'   all three censoring levels. PB-only uses two at QRP none and all three at
#'   positive QRP, since those conditions are all biased. Both universes have
#'   the same groups. `"none"` retains full individual conditions.
#'   Each condition's mean tau averages valid estimates among its 1,000
#'   replications; `n_tau` records the count. `mean_tau2` is the average estimated
#'   variance, not the square of `mean_tau`.
#'
#' @details If K was simulated, only that K is searched. Otherwise the nearest K
#'   below and above are searched separately. At each K, the closest mean tau
#'   below and above the target is returned. Exact tau matches take priority;
#'   all ties within numerical tolerance (1e-12 on the tau scale) are retained.
#'   Outside the simulated range, only the available boundary is returned.
#'
#'   Rejection rates exclude missing tests. Applicability is valid replications
#'   divided by total replications. A condition is biased when censoring or QRP
#'   is non-zero. TPR is available for biased conditions, FPR for unbiased
#'   conditions; Balanced Accuracy requires both. Missing metrics remain NA.
#'   Pooled rates give equal weight to component-condition rates.
#'
#'   The returned summary gives equal weight to each selected match. A method's
#'   mean Balanced Accuracy or rejection rate is available only when that metric
#'   is finite for every match. Missing values are never replaced with zero.
#'   Summary rows are sorted by mean Balanced Accuracy, with unavailable means
#'   last. Low applicability is reported separately and retained in the output.
#'
#' @return An object of class `pbias_closest_result` with `matches` (condition
#'   settings), `performance` (numeric results for all 21 methods per match),
#'   `summary` (one row per method, averaged across all matches),
#'   `candidates` (the eligible conditions), `notes` (boundary messages), and
#'   `meta` (inputs, aggregation and simulation-data provenance). Rates are
#'   proportions, not percentages. Printing shows the matches on
#'   both heterogeneity scales and the mean Balanced Accuracy table. Defaults
#'   are Extended, all QRP levels, all conditions and pooled censoring levels.
#' @export
#' @examples
#' near <- pbias_closest(k = 65, tau = 0.2)
#' near
#' head(near$performance)
#' # If your estimate is tau-squared = 0.04:
#' pbias_closest(k = 60, tau = sqrt(0.04), universe = "all")
pbias_closest <- function(k, tau, design = c("extended", "original"),
                          qrp = c("all", "without", "with"),
                          universe = c("all", "pb_only"),
                          pooling = c("publication_bias", "none")) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 3 || k != floor(k)) {
    stop("k must be one finite whole number of studies, at least 3.", call. = FALSE)
  }
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0) {
    stop("tau must be one finite, non-negative estimated standard deviation; use sqrt(tau2) for a variance.", call. = FALSE)
  }
  design <- match.arg(design)
  qrp <- match.arg(qrp)
  universe <- match.arg(universe)
  pooling <- match.arg(pooling)
  data <- readRDS(.pbias_extdata_file("closest_conditions.rds"))
  if (pooling == "none") {
    candidates <- data$conditions
    candidates$profile_id <- paste0("condition-", candidates$Environment)
    candidates$n_conditions <- 1L
    candidates$n_unbiased_conditions <- as.integer(!candidates$biased)
    candidates$n_biased_conditions <- as.integer(candidates$biased)
    if (universe == "pb_only") candidates <- candidates[candidates$biased, , drop = FALSE]
    performance <- data$condition_performance
    performance$profile_id <- paste0("condition-", performance$Environment)
    index <- match(performance$Environment, data$conditions$Environment)
    performance$n_total <- data$conditions$n_total[index]
    biased <- data$conditions$biased[index]
    performance$rejection_rate <- performance$rate
    performance$applicability <- performance$n_valid / performance$n_total
    performance$tpr <- ifelse(biased, performance$rate, NA_real_)
    performance$fpr <- ifelse(!biased, performance$rate, NA_real_)
    performance$balanced_accuracy <- NA_real_
  } else {
    candidates <- data$profiles[data$profiles$universe == universe, , drop = FALSE]
    candidates$Censor <- ifelse(universe == "pb_only" & candidates$QRP == "none",
                                 "med/high pooled", "none/med/high pooled")
    candidates$Environment <- NA_integer_
    performance <- data$performance
  }
  if (design == "original") candidates <- candidates[candidates$K >= 10 & candidates$K <= 100, , drop = FALSE]
  if (qrp == "without") candidates <- candidates[candidates$QRP == "none", , drop = FALSE]
  if (qrp == "with") candidates <- candidates[candidates$QRP != "none", , drop = FALSE]
  candidates <- candidates[order(candidates$K, candidates$mean_tau, candidates$profile_id), , drop = FALSE]
  rownames(candidates) <- NULL
  matches <- .pbias_bracket(candidates, k, tau)
  matches$match_id <- paste0("M", seq_len(nrow(matches)))
  matches$match_label <- paste0(matches$match_id, ": K = ", matches$K,
                                ", mean tau = ", formatC(matches$mean_tau, format = "f", digits = 5))
  perf <- performance[performance$profile_id %in% matches$profile_id, , drop = FALSE]
  perf$match_id <- matches$match_id[match(perf$profile_id, matches$profile_id)]
  perf$status <- ifelse(perf$n_valid == 0, "Unavailable",
                        ifelse(perf$applicability < 0.9, "Below 90% applicability", "Meets 90% threshold"))
  methods <- c("BEGG", paste0("EGGER", 1:5), "EK", "SKEW1", "SKEW2", "PSM3", "PSM4",
               "AK1", "AK2", "PUNIF", "TES", "PSST", "TESS", "CALI05", "CALI10", "CALI15", "CALI20")
  perf <- perf[order(match(perf$match_id, matches$match_id), match(perf$method, methods)),
               c("match_id", "profile_id", "method", "rejection_rate", "applicability", "n_valid", "n_total",
                 "tpr", "fpr", "balanced_accuracy", "status"), drop = FALSE]
  rownames(perf) <- NULL
  summary <- .pbias_closest_summary(perf, matches$match_id)
  notes <- character()
  ks <- sort(unique(candidates$K))
  if (k < min(ks) || k > max(ks)) notes <- c(notes, paste0(
    "K is outside the simulated range (", min(ks), "-", max(ks), "). Only the nearest available K is shown."))
  for (value_k in unique(matches$K)) {
    values <- candidates$mean_tau[candidates$K == value_k]
    if (tau < min(values) || tau > max(values)) notes <- c(notes, paste0(
      "At K = ", value_k, ", tau is outside the available range; only the nearest boundary condition is shown."))
  }
  structure(list(matches = matches, performance = perf, summary = summary, candidates = candidates, notes = notes,
                 meta = c(list(k = k, tau = tau, design = design, qrp = qrp, universe = universe,
                               pooling = pooling, package_version = "0.3.1",
                               summary_weighting = "Equal weight per match; all matches must have a valid metric"), data$metadata)),
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
  cat("Closest simulation conditions\n")
  cat(sprintf("K = %s; estimated tau = %.6g; tau-squared = %.6g\n",
              x$meta$k, x$meta$tau, x$meta$tau^2))
  cat(sprintf("%s design | %s QRP | %s universe | %s pooling\n",
              x$meta$design, x$meta$qrp, x$meta$universe, x$meta$pooling))
  m <- x$matches
  print(data.frame(Match = m$match_id, K = m$K, Mean_tau = round(m$mean_tau, 6),
                   Mean_tau2 = round(m$mean_tau2, 6),
                   Side = sub("Tau ", "", sub("Exact tau", "Exact", m$heterogeneity_match))), row.names = FALSE)
  cat("Full condition settings: $matches. All 21 methods per match: $performance.\n\n")
  s <- x$summary
  pct <- function(z) ifelse(is.finite(z), sprintf("%.1f%%", 100 * z), "n/a")
  print(data.frame(Method = s$method, Mean_BA = pct(s$mean_balanced_accuracy),
                   BA_matches = paste0(s$ba_available, "/", s$n_matches),
                   Mean_applicability = pct(s$mean_applicability),
                   Low_applicability = s$low_applicability_matches), row.names = FALSE)
  cat("Equal weight per match; an average requires valid values in every match.\n")
  cat("Low applicability counts matches below 90%; inspect these before comparing methods.\n")
  if (length(x$notes)) cat(paste(x$notes, collapse = "\n"), "\n")
  invisible(x)
}

.pbias_closest_summary <- function(performance, match_ids) {
  methods <- unique(performance$method)
  n <- length(match_ids)
  complete_mean <- function(z) if (length(z) == n && all(is.finite(z))) mean(z) else NA_real_
  out <- do.call(rbind, lapply(methods, function(method) {
    p <- performance[performance$method == method, , drop = FALSE]
    p <- p[match(match_ids, p$match_id), , drop = FALSE]
    data.frame(method = method, mean_balanced_accuracy = complete_mean(p$balanced_accuracy),
               mean_rejection_rate = complete_mean(p$rejection_rate),
               mean_applicability = complete_mean(p$applicability),
               min_applicability = if (all(is.finite(p$applicability))) min(p$applicability) else NA_real_,
               n_matches = n, ba_available = sum(is.finite(p$balanced_accuracy)),
               rejection_available = sum(is.finite(p$rejection_rate)),
               low_applicability_matches = sum(!is.finite(p$applicability) | p$applicability < 0.9))
  }))
  out <- out[order(-out$mean_balanced_accuracy, match(out$method, methods), na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}
