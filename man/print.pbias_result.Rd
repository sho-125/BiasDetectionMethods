\name{print.pbias_result}
\alias{print.pbias_result}
\title{Print a pbiasr result}
\description{
Prints a styled summary with the null hypothesis, paper case, and the three-column publication-bias table ranked by Balanced Accuracy: Estimator, P-value, and Balanced Accuracy. The bottom note explains the columns, marks unavailable case-specific estimators as not recommended, and points users to the Shiny app.
}
\usage{
\method{print}{pbias_result}(x, digits = 3, ...)
}
\arguments{
  \item{x}{An object returned by \code{pbias_table()}.}
  \item{digits}{Number of digits to print.}
  \item{...}{Unused.}
}
\value{
Returns \code{x} invisibly.
}
