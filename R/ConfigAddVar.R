ConfigAddVar <- function(configuration, var_type, var_entry) {
  if (var_type == '2d') {
    if (!(var_entry %in% configuration$two_d_vars)) {
      configuration$two_d_vars <- c(configuration$two_d_vars, var_entry)
    }
  } else if (var_type == 'mean') {
    if (!(var_entry %in% configuration$global_mean_vars)) {
      configuration$global_mean_vars <- c(configuration$global_mean_vars, var_entry)
    }
  } else {
    stop("'var_type' must be either '2d' or 'mean'")
  }

  configuration
}
