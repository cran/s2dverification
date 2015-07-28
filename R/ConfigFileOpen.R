ConfigFileOpen <- function(file_path, silent = FALSE) {
  if (!silent) {
    cat(paste("* Reading configuration file:", file_path, "\n"))
  }
  # Read the data from the configuration file.
  ## Remove comments, tabulations, spaces, empty lines, ...
  all_lines <- readLines(file_path)
  all_lines <- gsub("\t", "", all_lines)
  all_lines <- gsub(" ", "", all_lines)
  all_lines <- all_lines[-grep("^#", all_lines)]
  all_lines <- all_lines[-grep("^$", all_lines)]
  ## Detect key lines
  key_positions <- grep("^!!", all_lines)

  ## Check that the format of the configuration file is right.
  if (length(key_positions) != 12) {
    stop('Error: The configuration file is corrupted or outdated: the key lines do not match the expected pattern.')
  }

  ## Start parsing the configuration.
  ## Parse the 2-dimensional (atmospheric) and global mean (oceanic) variable 
  ## entries and keep them in two vectors
  if (key_positions[1] + 1 < key_positions[2]) {
    two_d_var_lines <- all_lines[(key_positions[1] + 1):(key_positions[2] - 1)]
    two_d_vars <- unlist(strsplit(two_d_var_lines, ","))
  } else {
    two_d_vars <- c()
  }
  if (key_positions[2] + 1 < key_positions[3]) {
    global_mean_var_lines <- all_lines[(key_positions[2] + 1):(key_positions[3] - 1)]
    global_mean_vars <- unlist(strsplit(gsub(" ", "", global_mean_var_lines), ","))
  } else {
    global_mean_vars <- c()
  }

  # The variables that are used in the configuration filed are kept in 
  # 'definitions', an associative array (key-value array or dictionary).
  definitions <- list()
  ## Parse the variables definitions in the whole configuration file
  if (key_positions[3] + 1 < key_positions[4]) {
    exp_table_definitions <- all_lines[(key_positions[3] + 1):(key_positions[4] - 1)]
  } else {
    exp_table_definitions <- c()
  }
  if (key_positions[5] + 1 < key_positions[6]) {
    file_per_member_exp_table_definitions <- all_lines[(key_positions[5] + 1):(key_positions[6] - 1)]
  } else {
    file_per_member_exp_table_definitions <- c()
  }
  if (key_positions[7] + 1 < key_positions[8]) {
    obs_table_definitions <- all_lines[(key_positions[7] + 1):(key_positions[8] - 1)]
  } else {
    obs_table_definitions <- c()
  }
  if (key_positions[9] + 1 < key_positions[10]) {
    file_per_member_obs_table_definitions <- all_lines[(key_positions[9] + 1):(key_positions[10] - 1)]
  } else {
    file_per_member_obs_table_definitions <- c()
  }
  if (key_positions[11] + 1 < key_positions[12]) {
    file_per_dataset_obs_table_definitions <- all_lines[(key_positions[11] + 1):(key_positions[12] - 1)]
  } else {
    file_per_dataset_obs_table_definitions <- c()
  }
  all_definitions <- c(exp_table_definitions, file_per_member_exp_table_definitions, obs_table_definitions)
  all_definitions <- c(all_definitions, file_per_member_obs_table_definitions, file_per_dataset_obs_table_definitions)
  if (length(grep("=", all_definitions)) == length(all_definitions)) {
    for (definition in all_definitions) {
      if (length(which(strsplit(definition, "")[[1]] == "=")) == 1) {
        var_name <- strsplit(definition, "=")[[1]][1]
        tmp_value <- strsplit(definition, "=")[[1]][2]
        var_value <- ifelse(is.na(tmp_value), "", tmp_value)
        if ((length(which(strsplit(var_value, "")[[1]] == "$")) %% 2) == 0) {
          definitions[[var_name]] <- var_value
        } else {
          stop('Error: The configuration file is corrupted: there are incorrect variable definition lines in the definition zones. A closing "$" symbol may be missing.')
        }
      } else {
        stop('Error: The configuration file is corrupted: there are incorrect definition lines in the definition zones.')
      }
    }
  } else {
    stop('Error: The configuration file is corrupted: there are malformed definition lines in the definition zones.')
  }
  mandatory_definitions <- c("DEFAULT_EXP_MAIN_PATH", "DEFAULT_EXP_FILE_PATH", "DEFAULT_GRID", 
                             "DEFAULT_NC_VAR_NAME", "DEFAULT_SUFFIX", "DEFAULT_VAR_MIN", 
                             "DEFAULT_VAR_MAX", "EXP_FULL_PATH", "DEFAULT_OBS_MAIN_PATH", 
                             "DEFAULT_OBS_FILE_PATH", "OBS_FULL_PATH", "DEFAULT_DIM_NAME_LONGITUDES",
                             "DEFAULT_DIM_NAME_LATITUDES", "DEFAULT_DIM_NAME_MEMBERS")
  if (any(!(mandatory_definitions %in% names(definitions)))) {
    cat("Error: The following variables must be defined in the configuration file. You can add them with ConfigFileOpen(), ConfigEditDefinition() and ConfigFileSave() or by editing the configuration file by hand, as specified in ?ConfigFileOpen.")
    stop(paste(mandatory_definitions, collapse = ', '))
  }

  # Parse the entries in the tables
  ## These are the indices of the key positions in the vector of key positions
  tables_key_positions <- c(4, 6, 8, 10, 12)
  current_table <- 1
  for (table_key_position in tables_key_positions) {
    datasets <- list(c(), c(), c(), c())

    if (current_table < 3) {
      default_values <- c("$DEFAULT_EXP_MAIN_PATH$", "$DEFAULT_EXP_FILE_PATH$", "$DEFAULT_GRID$", "$DEFAULT_NC_VAR_NAME$", '$DEFAULT_SUFFIX$', '$DEFAULT_VAR_MIN$', '$DEFAULT_VAR_MAX$')
    } else {
      default_values <- c("$DEFAULT_OBS_MAIN_PATH$", "$DEFAULT_OBS_FILE_PATH$", "$DEFAULT_NC_VAR_NAME$", '$DEFAULT_SUFFIX$', '$DEFAULT_VAR_MIN$', '$DEFAULT_VAR_MAX$')
    }
    previous_values <- c(".*", ".*", default_values)
    table_lines <- c()
    table_end <- ifelse(table_key_position == max(tables_key_positions), length(all_lines), key_positions[table_key_position + 1] - 1)
    if ((key_positions[table_key_position] + 1) <= table_end) {
      table_lines <- all_lines[(key_positions[table_key_position] + 1):table_end]
      table_lines <- strsplit(table_lines, ",")
    }

    for (entry in table_lines) {
      if (entry[1] == '"') {
        entry[1] <- previous_values[1]
      }
      if ((length(entry) > 1)) {
        if (entry[2] == '"') {
          entry[2] <- previous_values[2]
        }
      } else {
        stop('Error: The variable column must be defined in all the entries in the tables in the configuration file.')
      }
      for (value_position in 1:length(default_values)) {
        if ((length(entry) > value_position + 1)) {
          if (entry[value_position + 2] == '"') {
            entry[value_position + 2] <- previous_values[value_position + 2]
          }
        } else {
          entry[value_position + 2] <- '*'
        }
      }
      if (entry[1] == '.*') {
        if (entry[2] == '.*') {
          datasets[[1]] <- c(datasets[[1]], list(entry))
        } else {
          datasets[[3]] <- c(datasets[[3]], list(entry))
        }
      } else {
        if (entry[2] == '.*') {
          datasets[[2]] <- c(datasets[[2]], list(entry))
        } else {
          datasets[[4]] <- c(datasets[[4]], list(entry))
        }
      }
      previous_values <- entry
    }

    if (current_table == 1) {
      file_per_startdate_exps <- datasets
    } else if (current_table == 2) {
      file_per_member_exps <- datasets
    } else if (current_table == 3) {
      file_per_month_obs <- datasets
    } else if (current_table == 4) {
      file_per_member_obs <- datasets
    } else if (current_table == 5) {
      file_per_dataset_obs <- datasets
    }

    current_table <- current_table + 1
  }
  
  if (!silent) {
    cat("* Config file read successfully.\n")
  }

  invisible(list(definitions = definitions, two_d_vars = two_d_vars, 
                 global_mean_vars = global_mean_vars, 
                 file_per_startdate_experiments = file_per_startdate_exps, 
                 file_per_member_experiments = file_per_member_exps, 
                 file_per_month_observations = file_per_month_obs,
                 file_per_member_observations = file_per_member_obs,
                 file_per_dataset_observations = file_per_dataset_obs))
}
