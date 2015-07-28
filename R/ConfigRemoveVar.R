ConfigRemoveVar <- function(configuration, var_type, var_entry) {
  if (var_type == '2d') {
    if (var_entry %in% configuration$two_d_vars) {
      configuration$two_d_vars <- configuration$two_d_vars[- match(var_entry, configuration$two_d_vars)]
    }
  } else if (var_type == 'mean') {
    if (var_entry %in% configuration$global_mean_vars) {
      configuration$global_mean_vars <- configuration$global_mean_vars[- match(var_entry, configuration$global_mean_vars)]
    }
  } else {
    stop("'var_type' must be either '2d' or 'mean'")
  }
  
  configuration
}
