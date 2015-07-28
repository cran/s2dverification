ConfigShowVars <- function(configuration) {
  cat("---------------------------\n")
  cat("| 2-dimensional variables |\n")
  cat("---------------------------\n")
  invisible(lapply(paste0(configuration$two_d_vars, "\n"), cat))
  cat("\n")
  cat("-------------------------\n")
  cat("| global-mean variables |\n")
  cat("-------------------------\n")
  invisible(lapply(paste0(configuration$global_mean_vars, "\n"), cat))
}
