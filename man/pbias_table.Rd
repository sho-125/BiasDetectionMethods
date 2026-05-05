\name{pbias_table}
\alias{pbias_table}
\title{Run publication-bias methods for a meta-analysis}
\description{
Runs the unified publication-bias workflow, classifies the meta-analysis into the paper's size-by-I-squared case, and returns a three-column publication-bias table sorted by the corresponding Balanced Accuracy ranking.
}
\usage{
pbias_table(
  yi,
  sei,
  obs = NULL,
  study_id = NULL,
  alpha = 0.05,
  sort_by_rank = TRUE
)
}
\arguments{
  \item{yi}{Numeric vector of observed effect sizes.}
  \item{sei}{Numeric vector of standard errors associated with \code{yi}.}
  \item{obs}{Optional numeric vector of study-level sample size proxies. Required for EGGER4 / EGGER5 in this implementation.}
  \item{study_id}{Optional study identifiers. Defaults to \code{1:length(yi)}.}
  \item{alpha}{Significance level used by the caliper windows and excess-significance tests.}
  \item{sort_by_rank}{Logical. If \code{TRUE}, sort the displayed publication-bias results by the Balanced Accuracy ranking for the detected paper case.}
}
\value{
An object of class \code{pbias_result} with elements \code{bias_summary}, \code{publication_table}, \code{case}, \code{interpretation}, and \code{meta}. \code{publication_table} has columns \code{Estimator}, \code{P-value}, and \code{Balanced Accuracy}. \code{bias_summary} has no row names, omits the internal code-method column, and labels the raw publication-bias p-value column as \code{p-value}. Estimated effects and model-parameter estimates are not reported.
}
\examples{
yi <- c(0.20, 0.35, 0.10, 0.42, 0.30, 0.28)
sei <- c(0.10, 0.12, 0.15, 0.11, 0.09, 0.13)
out <- pbias_table(yi, sei)
print(out)
out$case
out$publication_table
out$bias_summary
}
