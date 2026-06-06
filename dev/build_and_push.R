# ============================================================
# update_package_and_push.R
# Update R package, build .tar.gz, install locally,
# commit changes, and push to GitHub
# ============================================================

# ---------- 0. Check you are in the package folder ----------
cat("\nCurrent working directory:\n")
print(getwd())

cat("\nFiles in this folder:\n")
print(list.files())

if (!file.exists("DESCRIPTION")) {
  stop("DESCRIPTION file not found. You are not in the R package root folder.")
}

# ---------- 1. Install required packages if missing ----------
needed_packages <- c("devtools", "remotes")

for (pkg in needed_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# ---------- 2. Set Git identity ----------
system('git config --global user.name "sho-125"')
system('git config --global user.email "sanghyun.hong@canterbury.ac.nz"')

# Also set it for this repo specifically
system('git config user.name "sho-125"')
system('git config user.email "sanghyun.hong@canterbury.ac.nz"')

# ---------- 3. Check Git remote ----------
cat("\nGit remote:\n")
system("git remote -v")

# If this shows nothing, uncomment and run this line once:
# system("git remote add origin https://github.com/sho-125/BiasDetectionMethods.git")

# ---------- 4. Update .Rbuildignore ----------
# This prevents temporary/build files from being included inside the R package build.

if (!file.exists(".Rbuildignore")) {
  file.create(".Rbuildignore")
}

rbuildignore <- readLines(".Rbuildignore", warn = FALSE)

needed_ignore <- c(
  "^dist$",
  "^.*\\.Rcheck$",
  "^update_package_and_push\\.R$",
  "^\\.Rproj\\.user$",
  "^\\.Rhistory$",
  "^\\.RData$"
)

for (x in needed_ignore) {
  if (!any(rbuildignore == x)) {
    cat(x, file = ".Rbuildignore", append = TRUE, sep = "\n")
  }
}

# ---------- 5. Update documentation ----------
cat("\nUpdating documentation...\n")
devtools::document()

# If README.Rmd exists, rebuild README.md
if (file.exists("README.Rmd")) {
  cat("\nRebuilding README.md...\n")
  devtools::build_readme()
}

# ---------- 6. Check package ----------
cat("\nChecking package...\n")
devtools::check(error_on = "error")

# ---------- 7. Build source package .tar.gz ----------
cat("\nBuilding source package tar.gz...\n")

if (!dir.exists("dist")) {
  dir.create("dist")
}

tar_file <- devtools::build(path = "dist")

cat("\nBuilt source package:\n")
print(tar_file)

# ---------- 8. Install local built package ----------
cat("\nInstalling built package locally...\n")
install.packages(tar_file, repos = NULL, type = "source")

# ---------- 9. Confirm package version ----------
cat("\nInstalled pbiasr version:\n")
library(pbiasr)
print(packageVersion("pbiasr"))

# ---------- 10. Check Git status ----------
cat("\nGit status before commit:\n")
system("git status")

# ---------- 11. Add files to Git ----------
system("git add -A")

# Force-add the tar.gz file in dist/
# This is useful if your .gitignore ignores *.tar.gz
tar_file_for_git <- gsub("\\\\", "/", tar_file)
system(paste("git add -f", shQuote(tar_file_for_git)))

# ---------- 12. Commit changes ----------
commit_msg <- "Update pbiasr package"

changed_files <- system("git status --porcelain", intern = TRUE)

if (length(changed_files) == 0) {
  cat("\nNo changes to commit.\n")
} else {
  cat("\nCommitting changes...\n")
  system(paste0('git commit -m "', commit_msg, '"'))
}

# ---------- 13. Check branch ----------
cat("\nCurrent Git branch:\n")
system("git branch")

# ---------- 14. Push to GitHub ----------
cat("\nPushing to GitHub...\n")
system("git push -u origin main")

# If your branch is master instead of main, use this instead:
# system("git push -u origin master")

# ---------- 15. Final status ----------
cat("\nFinal Git status:\n")
system("git status")

cat("\nDone.\n")
cat("Package updated, documentation regenerated, checked, tar.gz built, installed locally, committed, and pushed to GitHub.\n")
cat("\nYour tar.gz file is here:\n")
print(tar_file)













# installing the package (need remotes package)
remotes::install_github(
  "sho-125/BiasDetectionMethods",
  dependencies = TRUE,
  upgrade = "always",
  force = TRUE
)


# Load package
library(pbiasr)
sample_data1 <- pbias_sample_data(1)
out1 <- pbias_table(
  yi = sample_data1$eff,
  sei = sample_data1$se,
  obs = sample_data1$obs
)
print(out1)
