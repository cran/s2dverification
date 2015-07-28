ConfigEditEntry <- function(configuration, dataset_type, store_format, position, dataset_name = NULL, var_name = NULL, main_path = NULL, file_path = NULL, nc_var_name = NULL, grid = NULL, suffix = NULL, varmin = NULL, varmax = NULL) {
  table_name <- paste0(gsub("-", "_", store_format), "_", dataset_type)
  
  all_entries <- length(unlist(configuration[[table_name]], recursive = FALSE))
  if (position < 1 || position > all_entries) {
    stop("'position' must be in the range [1, # of table entries]")
  }

  found <- FALSE
  level <- 1
  index_of_first <- 1
  while (!found && level < 5) {
    if (position <= (index_of_first + length(configuration[[table_name]][[level]]) - 1)) {
      found <- TRUE
    } else {
      index_of_first <- index_of_first + length(configuration[[table_name]][[level]])
      level <- level + 1
    }
  }
  position <- position - index_of_first + 1

  if (dataset_type == 'experiments') {
    edited_values <- c(1:9)[c(!is.null(dataset_name), !is.null(var_name), !is.null(main_path), !is.null(file_path), !is.null(grid), !is.null(nc_var_name), !is.null(suffix), !is.null(varmin), !is.null(varmax))]
    configuration[[table_name]][[level]][[position]][edited_values] <- c(dataset_name, var_name, main_path, file_path, grid, nc_var_name, suffix, varmin, varmax)
  } else if (dataset_type == 'observations') {
    edited_values <- c(1:8)[c(!is.null(dataset_name), !is.null(var_name), !is.null(main_path), !is.null(file_path), !is.null(nc_var_name), !is.null(suffix), !is.null(varmin), !is.null(varmax))]
    configuration[[table_name]][[level]][[position]][edited_values] <- c(dataset_name, var_name, main_path, file_path, nc_var_name, suffix, varmin, varmax)
  } else {
    stop("'dataset_type' must be one of 'experiments' or 'observations'")
  }

  configuration
}
