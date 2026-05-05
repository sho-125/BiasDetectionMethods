\name{pbias_sample_data}
\alias{pbias_sample_data}
\title{Load packaged sample data}
\description{
Loads one of the sample datasets from \code{inst/extdata}. Use \code{which = 1} for \code{SampleData1.RData}, \code{which = 2} for \code{SampleData2.RData}, or \code{which = "all"} to load both.
}
\usage{
pbias_sample_data(which = 1)
}
\arguments{
  \item{which}{Which sample dataset to load. Use \code{1}, \code{2}, \code{"SampleData1"}, \code{"SampleData2"}, or \code{"all"}.}
}
\value{
A sample-data object, or a named list of sample-data objects when \code{which = "all"} or when a file contains multiple objects.
}
\examples{
sample_data1 <- pbias_sample_data(1)
sample_data2 <- pbias_sample_data(2)
}
