stopifnot(file.exists("DESCRIPTION"))
stopifnot(file.exists("pbiasr.Rproj"))

# Install packages if needed
needed <- c("usethis", "devtools", "pkgbuild")
missing <- needed[!needed %in% rownames(installed.packages())]

if (length(missing) > 0) {
  install.packages(missing)
}

# Optional: bump patch version in DESCRIPTION
# Example: 0.1.0 -> 0.1.1
usethis::use_version("patch")

# Rebuild documentation if you use roxygen comments
devtools::document()

# Check package
devtools::check()

# Build source tar.gz file outside the package folder
build_dir <- file.path(dirname(getwd()), "build")
dir.create(build_dir, showWarnings = FALSE)

tarball <- pkgbuild::build(
  path = ".",
  dest_path = build_dir,
  binary = FALSE
)

cat("Built tarball:\n", tarball, "\n")

# Git setup
if (!dir.exists(".git")) {
  usethis::use_git()
}

system('git config user.name "sho-125"')
system('git config user.email "sanghyun.hong@canterbury.ac.nz"')

# Use main branch
system("git branch -M main")

# Set GitHub remote
remote_exists <- system(
  "git remote get-url origin",
  ignore.stdout = TRUE,
  ignore.stderr = TRUE
) == 0

if (remote_exists) {
  system("git remote set-url origin https://github.com/sho-125/BiasDetectionMethods.git")
} else {
  system("git remote add origin https://github.com/sho-125/BiasDetectionMethods.git")
}

# Add, commit, and push
system("git add .")
system("git status")
system('git commit -m "Show ldist in publication table"')
system("git push -u origin main")