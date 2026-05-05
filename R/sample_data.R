#' Load packaged sample data
#'
#' Loads one of the sample datasets from \code{inst/extdata}. Use
#' \code{which = 1} for \code{SampleData1.RData}, \code{which = 2} for
#' \code{SampleData2.RData}, or \code{which = "all"} to load both.
#'
#' @param which Which sample dataset to load. Use \code{1}, \code{2},
#'   \code{"SampleData1"}, \code{"SampleData2"}, or \code{"all"}.
#'
#' @return A sample-data object, or a named list of sample-data objects when
#'   \code{which = "all"} or when a file contains multiple objects.
#'
#' @examples
#' sample_data1 <- pbias_sample_data(1)
#' sample_data2 <- pbias_sample_data(2)
#'
#' @export
pbias_sample_data <- function(which = 1) {
  files <- c(SampleData1 = "SampleData1.RData", SampleData2 = "SampleData2.RData")

  if (identical(which, "all")) {
    return(stats::setNames(
      lapply(names(files), function(name) .load_sample_file(files[[name]])),
      names(files)
    ))
  }

  if (is.numeric(which) && length(which) == 1L && which %in% seq_along(files)) {
    return(.load_sample_file(files[[which]]))
  }

  if (is.character(which) && length(which) == 1L) {
    key <- sub("\\.RData$", "", which, ignore.case = TRUE)
    idx <- match(key, names(files))
    if (!is.na(idx)) {
      return(.load_sample_file(files[[idx]]))
    }
  }

  stop("`which` must be 1, 2, 'SampleData1', 'SampleData2', or 'all'.")
}

.load_sample_file <- function(file) {
  path <- system.file("extdata", file, package = "pbiasr")
  if (!nzchar(path) || !file.exists(path)) {
    local_path <- file.path("inst", "extdata", file)
    if (file.exists(local_path)) {
      path <- local_path
    }
  }
  if (!nzchar(path) || !file.exists(path)) {
    stop("Could not find `", file, "` in package extdata.")
  }

  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  if (length(object_names) == 1L) {
    return(get(object_names, envir = env))
  }

  stats::setNames(
    lapply(object_names, function(name) get(name, envir = env)),
    object_names
  )
}
