ConfigApplyMatchingEntries <- function(configuration, var, exp = NULL, obs = NULL, show_entries = FALSE, show_result = TRUE) {
  ## Function to tell if a regexpr() match is a complete match to a specified name
  isFullMatch <- function(x, name) {
    ifelse(x > 0 && attributes(x)$match.length == nchar(name), TRUE, FALSE)
  }

  ## This function calculates the maximum of a positive vector.
  positiveMax <- function(x) {
    max(c(0, x[!sapply(x, is.null) & !sapply(x, is.na)]))
  }

  var_entries_in_file_per_startdate_exps <- c()
  if (length(unlist(configuration$file_per_startdate_experiments, recursive = FALSE)) > 0) {
    var_entries_in_file_per_startdate_exps <- which(unlist(lapply(lapply(as.list(unlist(lapply(configuration$file_per_startdate_experiments, lapply, "[[", 2))), regexpr, var), isFullMatch, var) > 0))
  }
  var_entries_in_file_per_member_exps <- c()
  if (length(unlist(configuration$file_per_member_experiments, recursive = FALSE)) > 0) {
    var_entries_in_file_per_member_exps <- which(unlist(lapply(lapply(as.list(unlist(lapply(configuration$file_per_member_experiments, lapply, "[[", 2))), regexpr, var), isFullMatch, var) > 0))
  }
  var_entries_in_file_per_month_obs <- c()
  if (length(unlist(configuration$file_per_month_observations, recursive = FALSE)) > 0) {
    var_entries_in_file_per_month_obs <- which(unlist(lapply(lapply(as.list(unlist(lapply(configuration$file_per_month_observations, lapply, "[[", 2))), regexpr, var), isFullMatch, var) > 0))
  }
  var_entries_in_file_per_member_obs <- c()
  if (length(unlist(configuration$file_per_member_observations, recursive = FALSE)) > 0) {
    var_entries_in_file_per_member_obs <- which(unlist(lapply(lapply(as.list(unlist(lapply(configuration$file_per_member_observations, lapply, "[[", 2))), regexpr, var), isFullMatch, var) > 0))
  }
  var_entries_in_file_per_dataset_obs <- c()
  if (length(unlist(configuration$file_per_dataset_observations, recursive = FALSE)) > 0) {
    var_entries_in_file_per_dataset_obs <- which(unlist(lapply(lapply(as.list(unlist(lapply(configuration$file_per_dataset_observations, lapply, "[[", 2))), regexpr, var), isFullMatch, var) > 0))
  }

  exp_info <- list()
  is_file_per_member_exp <- rep(FALSE, length(exp))
  jmod <- 1
  for (mod in exp) {
    file_per_startdate_mod_var_matching_entries <- file_per_startdate_mod_var_matching_indices <- file_per_startdate_mod_var_matching_entries_levels <- c()
    file_per_member_mod_var_matching_entries <- file_per_member_mod_var_matching_indices <- file_per_member_mod_var_matching_entries_levels <- c()
    
    if (length(unlist(configuration$file_per_startdate_experiments, recursive = FALSE)) > 0) {
      file_per_startdate_mod_entries_in_exps <- which(unlist(lapply(lapply(unlist(lapply(configuration$file_per_startdate_experiments, lapply, "[[", 1), recursive = FALSE), regexpr, mod), isFullMatch, mod)))
      if (length(file_per_startdate_mod_entries_in_exps) > 0) {
        file_per_startdate_mod_var_matching_indices <- intersect(var_entries_in_file_per_startdate_exps, file_per_startdate_mod_entries_in_exps)
        file_per_startdate_mod_var_matching_entries <- unlist(configuration$file_per_startdate_experiments, recursive = FALSE)[file_per_startdate_mod_var_matching_indices]
        file_per_startdate_exps_levels <- lapply(as.list(1:4), f <- function(x) {x <- array(x, length(configuration$file_per_startdate_experiments[[x]]))})
        file_per_startdate_mod_var_matching_entries_levels <- unlist(file_per_startdate_exps_levels)[intersect(var_entries_in_file_per_startdate_exps, file_per_startdate_mod_entries_in_exps)]
      }
    }

    if (length(unlist(configuration$file_per_member_experiments, recursive = FALSE)) > 0) {
      file_per_member_mod_entries_in_exps <- which(unlist(lapply(lapply(unlist(lapply(configuration$file_per_member_experiments, lapply, "[[", 1), recursive = FALSE), regexpr, mod), isFullMatch, mod)))
      if (length(file_per_member_mod_entries_in_exps) > 0) {
        file_per_member_mod_var_matching_indices <- intersect(var_entries_in_file_per_member_exps, file_per_member_mod_entries_in_exps)
        file_per_member_mod_var_matching_entries <- unlist(configuration$file_per_member_experiments, recursive = FALSE)[file_per_member_mod_var_matching_indices]
        file_per_member_exps_levels <- lapply(as.list(1:4), f <- function(x) {x <- array(x, length(configuration$file_per_member_experiments[[x]]))})
        file_per_member_mod_var_matching_entries_levels <- unlist(file_per_member_exps_levels)[intersect(var_entries_in_file_per_member_exps, file_per_member_mod_entries_in_exps)]
      }
    }

    if (positiveMax(file_per_member_mod_var_matching_entries_levels) > positiveMax(file_per_startdate_mod_var_matching_entries_levels)) {
      is_file_per_member_exp[jmod] <- TRUE
      mod_var_matching_entries <- file_per_member_mod_var_matching_entries
    } else {
      mod_var_matching_entries <- file_per_startdate_mod_var_matching_entries
    }

    if (length(mod_var_matching_entries) == 0) {
      stop(paste('Error: There are no matching entries in the configuration file for the experiment', mod, 'and the variable', var, 
                 '. Please check the configuration file.)'))
    } else {
      if (show_entries) {
        header <- paste0("# Matching entries for experiment '", exp[jmod], "' and variable '", var, "' #\n")
        cat(paste0(paste(rep("#", nchar(header) - 1), collapse = ''), "\n"))
        cat(header)
        cat(paste0(paste(rep("#", nchar(header) - 1), collapse = ''), "\n"))
        cat(paste0("* (The entries that Load() will apply are the ones in ", 
          if (is_file_per_member_exp[jmod]) {
            "file per member"
          } else {
            "file per startdate"
          }, " experiments table, since contains most specific entries for the selected experiment and variable)\n"))
        ConfigShowTable(list(file_per_startdate_experiments = list(file_per_startdate_mod_var_matching_entries)), 'experiments', 'file-per-startdate', file_per_startdate_mod_var_matching_indices)
        ConfigShowTable(list(file_per_member_experiments = list(file_per_member_mod_var_matching_entries)), 'experiments', 'file-per-member', file_per_member_mod_var_matching_indices)
        cat("\n")
      }
      result <- .ConfigGetDatasetInfo(mod_var_matching_entries)
      if (show_result) {
        cat(paste0("The result of applying the matching entries to experiment name '", exp[jmod], "' and variable name '", var, "' is:\n"))
        configuration$definitions[["VAR_NAME"]] <- var
        configuration$definitions[["EXP_NAME"]] <- exp[jmod]
        fields <- c("MAIN_PATH: ", "FILE_PATH: ", "GRID: ", "NC_VAR_NAME: ", "SUFFIX: ", "VAR_MIN: ", "VAR_MAX: ")
        values <- lapply(result, lapply, function (x) .ConfigReplaceVariablesInString(x, configuration$definitions, TRUE))
        lapply(paste0(fields, unlist(values), "\n"), cat)
        cat("\n")
      }
      exp_info <- c(exp_info, list(result))
    }

    jmod <- jmod + 1
  }

  obs_info <- list()
  is_file_per_member_obs <- is_file_per_dataset_obs <- rep(FALSE, length(obs))
  jobs <- 1
  for (ref in obs) {
    file_per_month_ref_var_matching_entries <- file_per_month_ref_var_matching_indices <- file_per_month_ref_var_matching_entries_levels <- c()
    file_per_member_ref_var_matching_entries <- file_per_member_ref_var_matching_indices <- file_per_member_ref_var_matching_entries_levels <- c()
    file_per_dataset_ref_var_matching_entries <- file_per_dataset_ref_var_matching_indices <- file_per_dataset_ref_var_matching_entries_levels <- c()
    
    if (length(unlist(configuration$file_per_month_observations, recursive = FALSE)) > 0) {
      file_per_month_ref_entries_in_obs <- which(unlist(lapply(lapply(unlist(lapply(configuration$file_per_month_observations, lapply, "[[", 1), recursive = FALSE), regexpr, ref), isFullMatch, ref)))
      if (length(file_per_month_ref_entries_in_obs) > 0) {
        file_per_month_ref_var_matching_indices <- intersect(var_entries_in_file_per_month_obs, file_per_month_ref_entries_in_obs)
        file_per_month_ref_var_matching_entries <- unlist(configuration$file_per_month_observations, recursive = FALSE)[file_per_month_ref_var_matching_indices]
        file_per_month_obs_levels <- lapply(as.list(1:4), f <- function(x) {x <- array(x, length(configuration$file_per_month_observations[[x]]))})
        file_per_month_ref_var_matching_entries_levels <- unlist(file_per_month_obs_levels)[intersect(var_entries_in_file_per_month_obs, file_per_month_ref_entries_in_obs)]
      }
    }

    if (length(unlist(configuration$file_per_member_observations, recursive = FALSE)) > 0) {
      file_per_member_ref_entries_in_obs <- which(unlist(lapply(lapply(unlist(lapply(configuration$file_per_member_observations, lapply, "[[", 1), recursive = FALSE), regexpr, ref), isFullMatch, ref)))
      if (length(file_per_member_ref_entries_in_obs) > 0) {
        file_per_member_ref_var_matching_indices <- intersect(var_entries_in_file_per_member_obs, file_per_member_ref_entries_in_obs)
        file_per_member_ref_var_matching_entries <- unlist(configuration$file_per_member_observations, recursive = FALSE)[file_per_member_ref_var_matching_indices]
        file_per_member_obs_levels <- lapply(as.list(1:4), f <- function(x) {x <- array(x, length(configuration$file_per_member_observations[[x]]))})
        file_per_member_ref_var_matching_entries_levels <- unlist(file_per_member_obs_levels)[intersect(var_entries_in_file_per_member_obs, file_per_member_ref_entries_in_obs)]
      }
    }

    if (length(unlist(configuration$file_per_dataset_observations, recursive = FALSE)) > 0) {
      file_per_dataset_ref_entries_in_obs <- which(unlist(lapply(lapply(unlist(lapply(configuration$file_per_dataset_observations, lapply, "[[", 1), recursive = FALSE), regexpr, ref), isFullMatch, ref)))
      if (length(file_per_dataset_ref_entries_in_obs) > 0) {
        file_per_dataset_ref_var_matching_indices <- intersect(var_entries_in_file_per_dataset_obs, file_per_dataset_ref_entries_in_obs)
        file_per_dataset_ref_var_matching_entries <- unlist(configuration$file_per_dataset_observations, recursive = FALSE)[file_per_dataset_ref_var_matching_indices]
        file_per_dataset_obs_levels <- lapply(as.list(1:4), f <- function(x) {x <- array(x, length(configuration$file_per_dataset_observations[[x]]))})
        file_per_dataset_ref_var_matching_entries_levels <- unlist(file_per_dataset_obs_levels)[intersect(var_entries_in_file_per_dataset_obs, file_per_dataset_ref_entries_in_obs)]
      }
    }

    if (positiveMax(file_per_member_ref_var_matching_entries_levels) > 
        max(positiveMax(file_per_month_ref_var_matching_entries_levels), positiveMax(file_per_dataset_ref_var_matching_entries_levels))) {
      is_file_per_member_obs[jobs] <- TRUE
      ref_var_matching_entries <- file_per_member_ref_var_matching_entries
    } else if (positiveMax(file_per_dataset_ref_var_matching_entries_levels) > 
        max(positiveMax(file_per_month_ref_var_matching_entries_levels), positiveMax(file_per_member_ref_var_matching_entries_levels))) {
      is_file_per_dataset_obs[jobs] <- TRUE
      ref_var_matching_entries <- file_per_dataset_ref_var_matching_entries
    } else {
      ref_var_matching_entries <- file_per_month_ref_var_matching_entries
    }

    if (length(ref_var_matching_entries) == 0) {
      stop(paste('Error: There are no matching entries in the configuration file for the observation', ref, 'and the variable', var, 
                 '. Please check the configuration file.)'))
    } else {
      if (show_entries) {
        header <- paste0("# Matching entries for observation '", obs[jobs], "' and variable '", var, "' #\n")
        cat(paste0(paste(rep("#", nchar(header) - 1), collapse = ''), "\n"))
        cat(header)
        cat(paste0(paste(rep("#", nchar(header) - 1), collapse = ''), "\n"))
        cat(paste0("* (The entries that Load() will apply are the ones in ", 
          if (is_file_per_member_obs[jobs]) {
            "file per member"
          } else if (is_file_per_member_obs[jobs]) {
            "file per dataset"
          } else {
            "file per month"
          }, " observations table, since contains most specific entries for the selected observation and variable)\n"))
        ConfigShowTable(list(file_per_month_observations = list(file_per_month_ref_var_matching_entries)), 'observations', 'file-per-month', file_per_month_ref_var_matching_indices)
        ConfigShowTable(list(file_per_member_observations = list(file_per_member_ref_var_matching_entries)), 'observations', 'file-per-member', file_per_member_ref_var_matching_indices)
        ConfigShowTable(list(file_per_dataset_observations = list(file_per_dataset_ref_var_matching_entries)), 'observations', 'file-per-dataset', file_per_dataset_ref_var_matching_indices)
        cat("\n")
      }
      result <- .ConfigGetDatasetInfo(ref_var_matching_entries)
      if (show_result) {
        cat(paste0("The result of applying the matching entries to observation name '", obs[jobs], "' and variable name '", var, "' is:\n"))
        configuration$definitions[['VAR_NAME']] <- var
        configuration$definitions[["OBS_NAME"]] <- obs[jobs]
        fields <- c("MAIN_PATH: ", "FILE_PATH: ", "NC_VAR_NAME: ", "SUFFIX: ", "VAR_MIN: ", "VAR_MAX: ")
        values <- lapply(result, lapply, function (x) .ConfigReplaceVariablesInString(x, configuration$definitions, TRUE))
        lapply(paste0(fields, unlist(values), "\n"), cat)
        cat("\n")
      }
      obs_info <- c(obs_info, list(result))
    }

    jobs <- jobs + 1
  }

  invisible(list(exp_info = exp_info, is_file_per_member_exp = is_file_per_member_exp, obs_info = obs_info, is_file_per_member_obs = is_file_per_member_obs, is_file_per_dataset_obs = is_file_per_dataset_obs))
}
