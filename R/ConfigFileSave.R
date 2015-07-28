ConfigFileSave <- function(configuration, file_path, confirm = TRUE) {
  continue <- TRUE
  if (file.exists(file_path)) {
    if (confirm) {
      while (continue != 'y' && continue != 'n') {
        continue <- readline(paste0("WARNING: The configuration file '", file_path, "' already exists. It will be replaced. Continue? (y/n)\n"))
      }
      continue <- ifelse(continue == 'y', TRUE, FALSE)
    }
  }
  if (continue) {
    file_conn <- file(file_path)
    file_text <- c(
"# s2dverification configuration file",
"#",
"# Check ?ConfigFileOpen after loading s2dverification for detailed ",
"# documentation on this configuration file.",
"",
"#########################",
"!!2-dimensional variables",
"#########################"
                  )
    file_text <- c(file_text, configuration$two_d_vars)
    file_text <- c(file_text,
"",
"#######################",
"!!global mean variables",
"#######################"
                  )
    file_text <- c(file_text, configuration$global_mean_vars)
   
    table_names <- c("file_per_startdate_experiments", "file_per_member_experiments", "file_per_month_observations", "file_per_member_observations", "file_per_dataset_observations")
    for (table_name in table_names) {
      if (table_name %in% c("file_per_startdate_experiments", "file_per_member_experiments")) {
        dataset_type <- 'exp'
      } else {
        dataset_type <- 'obs'
      }
      file_text <- c(file_text,
"",
  paste(rep("#", nchar(table_name) + 14), collapse = ''),
  paste0("!!", gsub("_", " ", table_name), " definitions"),
  paste(rep("#", nchar(table_name) + 14), collapse = '')
                  )
      if (table_name == "file_per_startdate_experiments") {
        defaults <- configuration$definitions[grep("^DEFAULT_", names(configuration$definitions))]
        definitions <- configuration$definitions[-grep("^DEFAULT_", names(configuration$definitions))]
        file_text <- c(file_text, as.vector(paste(names(defaults), unlist(defaults), sep = " = ")))
        file_text <- c(file_text, as.vector(paste(names(definitions), unlist(definitions), sep = " = ")))
      }
      file_text <- c(file_text,
"",
  paste(rep("#", nchar(table_name) + 11), collapse = ''),
  paste0("!!table of ", gsub("_", " ", table_name)),
  paste(rep("#", nchar(table_name) + 11), collapse = ''),
  paste0("#", dataset_type, "_name, var_name[, ", dataset_type, "_main_path[, ", dataset_type, "_file_path[, ", ifelse(dataset_type == 'exp', "grid[, ", ""), "nc_var_name[, suffix[, var_min[, var_max]]]", ifelse(dataset_type == 'exp', "]", ""), "]]]")
                  )
      # Some heavy entry processing still to do here, to put asterisks, empty spaces, double quotes, and reduce options
      file_text <- c(file_text, unlist(lapply(configuration[[table_name]], function (x) lapply(x, function (y) paste(unlist(y), collapse = ", ")))))
    }
  
    writeLines(file_text, file_conn)
    close(file_conn)
  }

  invisible(continue)
}
