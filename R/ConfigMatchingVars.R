ConfigMatchingVars <- function(configuration, var, show = TRUE) {
  ## Function to tell if a regexpr() match is a complete match to a specified name
  isFullMatch <- function(x, name) {
    ifelse(x > 0 && attributes(x)$match.length == nchar(name), TRUE, FALSE)
  }

  if (length(configuration$two_d_vars) > 0) {
    entries_in_2d_vars <- which(unlist(lapply(lapply(as.list(configuration$two_d_vars), regexpr, var), isFullMatch, var)))
    has_entries_in_2d_vars <- length(entries_in_2d_vars) > 0
  } else {
    has_entries_in_2d_vars <- FALSE
  }
  if (length(configuration$global_mean_vars) > 0) {
    entries_in_global_mean_vars <- which(unlist(lapply(lapply(as.list(configuration$global_mean_vars), regexpr, var), isFullMatch, var)))
    has_entries_in_global_mean_vars <- length(entries_in_global_mean_vars) > 0
  } else {
    has_entries_in_global_mean_vars <- FALSE
  }
  if (has_entries_in_2d_vars) {
    is_2d_var <- TRUE
    if (show) {
      cat(paste("* In list of 2-dimensional variables:", configuration$two_d_vars[(entries_in_2d_vars)[1]], "\n"))
    }
  } else if (has_entries_in_global_mean_vars) {
    is_2d_var <- FALSE
    if (show) {
      cat(paste("* In list of global mean variables:", configuration$global_mean_vars[(entries_in_global_mean_vars)[1]], "\n"))
    }
  } else {
    stop('Error: The variable specified is not available. Please check the supported variables in the configuration file.')
  }

  invisible(is_2d_var)
}
