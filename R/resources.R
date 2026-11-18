# A private source root allows the reviewed package code to be bundled with the
# app without requiring a public GitHub release or changing the system library.
.pbias_source_dir <- NULL

.pbias_extdata_file <- function(file) {
  if (!is.null(.pbias_source_dir)) {
    return(file.path(.pbias_source_dir, "inst", "extdata", file))
  }
  path <- system.file("extdata", file, package = "pbiasr")
  if (!nzchar(path) || !file.exists(path)) {
    local_path <- file.path("inst", "extdata", file)
    if (file.exists(local_path)) path <- local_path
  }
  path
}
