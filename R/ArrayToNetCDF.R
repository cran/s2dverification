ArrayToNetCDF <- function(arrays, file_path) {
  # Check parameter arrays.
  if (is.array(arrays)) {
    arrays <- list(arrays)
  }
  if (any(!sapply(arrays, function(x) is.array(x) && (is.numeric(x) || is.logical(x))))) {
    stop("The parameter 'arrays' must be one or a list of numeric or logical arrays.")
  }
  # Check parameter file_path.
  if (!is.character(file_path)) {
    stop("Parameter 'file_path' must be a character string.")
  }

  defined_dims <- list()
  defined_vars <- list()
  global_attrs <- list()
  var_dim <- NULL
  for (i in 1:length(arrays)) {
    array_attrs <- attributes(arrays[[i]])
    if ('variables' %in% names(array_attrs)) {
      vars_info <- array_attrs[['variables']]
      array_attrs <- array_attrs[-which(names(array_attrs) == 'variables')]
    } else {
      vars_info <- NULL
    }
    global_attrs[names(array_attrs)] <- array_attrs
    var_dim <- which(names(dim(arrays[[i]])) %in% c('var', 'variable'))
    if (length(var_dim) > 0) {
      var_dim <- var_dim[1]
      num_vars <- dim(arrays[[i]])[var_dim]
    } else {
      var_dim <- NULL
      num_vars <- 1
    }
    # Defining ncdf4 variable objects
    for (j in 1:num_vars) {
      var_info <- vars_info[[j]]
      if (length(var_info) == 0) {
        var_info <- list()
      }
      dim_names <- names(dim(arrays[[i]]))
      if (!is.null(dim_names)) {
        if (any(is.na(dim_names) || (sapply(dim_names, nchar) == 0))) {
          stop("The provided arrays must have all named dimensions or ",
               "all unnamed dimensions.")
        }
      }
      provided_dims <- sapply(var_info$dim, '[[', 'name')
      var_built_dims <- NULL
      for (k in 1:length(dim(arrays[[i]]))) {
        if (!identical(k, var_dim)) {
          final_dim_position <- k - ifelse(!is.null(var_dim) && k > var_dim, 1, 0)
          dim_name <- dim_names[k]
          if (!is.null(dim_name) && (dim_name %in% provided_dims)) {
            dim_info <- var_info$dim[[which(provided_dims == dim_name)]]
          } else {
            dim_info <- list()
          }
          if (!('name' %in% names(dim_info))) {
            if (!is.null(dim_name)) {
              dim_info[['name']] <- dim_name
            } else {
              dim_info[['name']] <- paste0('dim', final_dim_position)
            }
          } else {
            if (!is.character(dim_info[['name']])) {
              stop("The provided 'name' for the ", k, "th dimension in the ", i, "th array must be a character string.")
            }
            dim_info[['name']] <- dim_info[['name']][1]
          }
          if (!('len' %in% names(dim_info))) {
            dim_info[['len']] <- unname(dim(arrays[[i]])[k])
          } else {
            if (!is.numeric(dim_info[['len']])) {
              stop("The provided 'len' for the ", k, "th dimension in the ", i, "th array must be a numeric value.")
            }
            dim_info[['len']] <- round(dim_info[['len']][1])
            if (dim_info[['len']] != dim(arrays[[i]])[k]) {
              stop("The provided 'len' for the ", k, "th dimension in the ", i, "th array does not match the actual length of the provided array.")
            }
          }
          if (!('unlim' %in% names(dim_info))) {
            dim_info[['unlim']] <- ifelse(dim_info[['name']] == 'time', TRUE, FALSE)
          } else {
            if (!is.logical(dim_info[['unlim']])) {
              stop("The provided 'unlim' for the ", k, "th dimension in the ", i, "th array must be a logical value.")
            }
            dim_info[['unlim']] <- dim_info[['unlim']][1]
          }
          if (!('units' %in% names(dim_info))) {
            dim_info[['units']] <- ''
          } else {
            if (!is.character(dim_info[['units']])) {
              stop("The provided 'units' for the ", k, "th dimension in the ", i, "th array must be a character string.")
            }
            dim_info[['units']] <- dim_info[['units']][1]
          }
          if (!('vals' %in% names(dim_info))) {
            dim_info[['vals']] <- 1:dim_info[['len']]
          } else {
            if (!is.numeric(dim_info[['vals']])) {
              stop("The provided 'vals' for the ", k, "th dimension in the ", i, "th array must be a numeric vector.")
            }
            if (dim_info[['units']] == '') {
              dim_info[['vals']] <- as.integer(dim_info[['vals']])
            }
            if (length(dim_info[['vals']]) != dim_info[['len']]) {
              stop("The length of the provided 'vals' for the ", k, "th dimension in the ", i, "th array does not match the length of the provided array.")
            }
          }
          if (!('create_dimvar' %in% names(dim_info))) {
            if (dim_info[['units']] == '') {
              dim_info[['create_dimvar']] <- FALSE
            } else {
              dim_info[['create_dimvar']] <- TRUE
            }
          } else {
            if (!is.logical(dim_info[['create_dimvar']])) {
              stop("The provided 'create_dimvar' for the ", k, "th dimension in the ", i, "th array must be a logical value.")
            }
            dim_info[['create_dimvar']] <- dim_info[['create_dimvar']][1]
            if (dim_info[['units']] != '' && !dim_info[['create_dimvar']]) {
              stop("Provided 'units' for the ", k, "th dimension in the ", i, "th array but 'create_dimvar' set to FALSE.")
            }
          }
          if (!('calendar' %in% names(dim_info))) {
            dim_info[['calendar']] <- NA
          } else {
            if (!is.character(dim_info[['calendar']])) {
              stop("The provided 'calendar' for the ", k, "th dimension in the ", i, "th array must be a character string.")
            }
            dim_info[['calendar']] <- dim_info[['calendar']][1]
          }
          if (!('longname' %in% names(dim_info))) {
            dim_info[['longname']] <- dim_info[['name']]
          } else {
            if (!is.character(dim_info[['longname']])) {
              stop("The provided 'longname' for the ", k, "th dimension in the ", i, "th array must be a character string.")
            }
            dim_info[['longname']] <- dim_info[['longname']][1]
          }
          if (dim_info[['name']] %in% names(defined_dims)) {
            items_to_check <- c('name', 'len', 'unlim', 'units', 'vals',
                                'create_dimvar', 'longname')
            if (!identical(dim_info[items_to_check], 
                           defined_dims[[dim_info[['name']]]][items_to_check]) ||
                !(identical(dim_info[['calendar']], defined_dims[[dim_info[['name']]]][['calendar']]) || 
                  (is.na(dim_info[['calendar']]) && is.null(defined_dims[[dim_info[['name']]]][['calendar']])))) {
              stop("The dimension '", dim_info[['name']], "' is ",
                   "defined or used more than once in the provided ",
                   "data but the dimension specifications do not ",
                   "match.")
            }
          } else {
            new_dim <- list(ncdim_def(dim_info[['name']], dim_info[['units']], 
                                      dim_info[['vals']], dim_info[['unlim']], 
                                      dim_info[['create_dimvar']], 
                                      dim_info[['calendar']], 
                                      dim_info[['longname']]))
            names(new_dim) <- dim_info[['name']]
            defined_dims <- c(defined_dims, new_dim)
          }
          var_built_dims <- c(var_built_dims, dim_info[['name']])
        }
      }
      if (!('name' %in% names(var_info))) {
        var_name_from_md <- names(vars_info)[j]
        var_name_from_ar <- names(arrays)[i]
        if (is.character(var_name_from_md) && 
            !is.na(var_name_from_md) &&
            (nchar(var_name_from_md) > 0)) {
          var_name <- var_name_from_md
        } else if (is.character(var_name_from_ar) && 
                   !is.na(var_name_from_ar) &&
                   (nchar(var_name_from_ar) > 0)){
          var_name <- var_name_from_ar
        } else {
          var_name <- paste0('var', i, '_', j)
        }
        var_info[['name']] <- var_name
      } else {
        if (!is.character(var_info[['name']])) {
          stop("The provided 'name' for the ", j, "th variable in the ", i, "th array must be a character string.")
        }
        var_info[['name']] <- var_info[['name']][1]
      }
      if (!('units' %in% names(var_info))) {
        var_info[['units']] <- ''
      } else {
        if (!is.character(var_info[['units']])) {
          stop("The provided 'units' for the ", j, "th variable in the ", i, "th array must be a character string.")
        }
        var_info[['units']] <- var_info[['units']][1]
      }
      if (!('missval' %in% names(var_info))) {
        var_info[['missval']] <- NULL
      } else {
        if (!is.numeric(var_info[['missval']])) {
          stop("The provided 'missval' for the ", j, "th variable in the ", i, "th array must be a numeric value.")
        }
        var_info[['missval']] <- var_info[['missval']][1]
      }
      if (!('longname' %in% names(var_info))) {
        var_info[['longname']] <- var_info[['name']]
      } else {
        if (!is.character(var_info[['longname']])) {
          stop("The provided 'longname' for the ", j, "th variable in the ", i, "th array must be a character string.")
        }
        var_info[['longname']] <- var_info[['longname']][1]
      }
      if (!('prec' %in% names(var_info))) {
        var_info[['prec']] <- 'float'
      } else {
        if (!is.character(var_info[['prec']])) {
          stop("The provided 'prec' for the ", j, "th variable in the ", i, "th array must be a character string.")
        }
        var_info[['prec']] <- var_info[['prec']][1]
      }
      new_var <- list(ncvar_def(var_info[['name']], var_info[['units']], 
                                defined_dims[var_built_dims], 
                                var_info[['missval']], 
                                var_info[['longname']], var_info[['prec']]))
      names(new_var) <- var_info[['name']]
      defined_vars <- c(defined_vars, new_var)
    }
  }
  ncdf_object <- nc_create(file_path, defined_vars)
  var_counter <- 1
  # Putting the data and extra attributes.
  for (i in 1:length(arrays)) {
    array_attrs <- attributes(arrays[[i]])
    if ('variables' %in% names(array_attrs)) {
      vars_info <- array_attrs[['variables']]
    } else {
      vars_info <- NULL
    }
    var_dim <- which(names(dim(arrays[[i]])) %in% c('var', 'variable'))
    if (length(var_dim) > 0) {
      var_dim <- var_dim[1]
      num_vars <- dim(arrays[[i]])[var_dim]
    } else {
      var_dim <- NULL
      num_vars <- 1
    }
    for (j in 1:num_vars) {
      var_info <- vars_info[[j]]
      if (length(var_info) == 0) {
        var_info <- list()
      }
      if (!('scaleFact' %in% names(var_info))) {
        scale_factor <- 1
      } else {
        if (!is.numeric(var_info[['scaleFact']])) {
          stop("The provided 'scaleFact' for the ", j, "th variable in the ", i, "th array must be a numeric value.")
        }
        scale_factor <- var_info[['scaleFact']][1]
      }
      if (!('addOffset' %in% names(var_info))) {
        add_offset <- 0
      } else {
        if (!is.numeric(var_info[['addOffset']])) {
          stop("The provided 'addOffset' for the ", j, "th variable in the ", i, "th array must be a numeric value.")
        }
        add_offset <- var_info[['addOffset']][1]
      }
      if (is.null(var_dim)) {
        if (scale_factor != 1 || add_offset != 0) {
          ncvar_put(ncdf_object, defined_vars[[var_counter]]$name, 
                    (arrays[[i]] - add_offset) / scale_factor, 
                    count = dim(arrays[[i]]))
        } else {
          ncvar_put(ncdf_object, defined_vars[[var_counter]]$name, 
                    arrays[[i]], 
                    count = dim(arrays[[i]]))
        }
      } else {
        if (scale_factor != 1 || add_offset != 0) {
          ncvar_put(ncdf_object, defined_vars[[var_counter]]$name, 
                    (Subset(arrays[[i]], var_dim, j, drop = 'selected') - add_offset) / scale_factor, 
                    count = dim(arrays[[i]])[-var_dim])
        } else {
          ncvar_put(ncdf_object, defined_vars[[var_counter]]$name, 
                    Subset(arrays[[i]], var_dim, j, drop = 'selected'), 
                    count = dim(arrays[[i]])[-var_dim])
        }
      }
      if (scale_factor != 1 || add_offset != 0) {
        ncatt_put(ncdf_object, defined_vars[[var_counter]]$name, 'scale_factor', scale_factor)
        ncatt_put(ncdf_object, defined_vars[[var_counter]]$name, 'add_offset', add_offset)
      }
      if ('coordinates' %in% names(var_info)) {
        if (!is.character(var_info[['coordinates']])) {
          stop("The attribute 'coordinates' must be a character string.")
        }
        if (!(all(strsplit(var_info[['coordinates']], ' ')[[1]] %in% sapply(defined_vars, '[[', 'name')))) {
          stop("All the dimensions appearing in 'coordinates' must point to defined variables.")
        }
        ncatt_put(ncdf_object, defined_vars[[var_counter]]$name, 'coordinates', var_info[['coordinates']])
      }
      var_counter <- var_counter + 1
    }
  }
  nc_close(ncdf_object)
  invisible(NULL)
}
