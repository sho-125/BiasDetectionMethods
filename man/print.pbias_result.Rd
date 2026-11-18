\name{print.pbias_result}
\alias{print.pbias_result}
\title{Print a pbiasr result}
\description{
Prints a styled summary with the null hypothesis, paper case, and the four-column publication-bias table ranked by the selected metric: Method, P-value, Balanced Accuracy, and Logit Distance. The bottom note explains the columns, marks unavailable case-specific estimators as Unranked*, and points users to the Shiny app.
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
