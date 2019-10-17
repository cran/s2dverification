#'Save multidimensional R arrays into NetCDF files
#'
#'This function takes as input one or a list of multidimensional R arrays and 
#'stores them in a NetCDF file, using the \code{ncdf4} package. The full path 
#'and name of the resulting file must be specified. Metadata can be attached 
#'to the arrays and propagated into the NetCDF file in 3 possible ways:
#'\itemize{
#'  \item{
#'  Via the list names if a list of arrays is provided: Each name in 
#'  the input list, corresponding to one multidimensional array, will be 
#'  interpreted as the name of the variable it contains.\cr
#'  E.g:\cr
#'  \code{ArrayToNetCDF(arrays = list(temperature = array(1:9, c(3, 3))),}\cr
#'  \code{                       file_path = 'example.nc')}
#'  }
#'  \item{
#'  Via the dimension names of each provided array: The dimension names 
#'  of each of the provided arrays will be interpreted as names for the 
#'  dimensions of the NetCDF files. Read further for special dimension names 
#'  that will trigger special behaviours, such as 'time' and 'var'.\cr
#'  E.g:\cr
#'  \code{temperature <- array(rnorm(100 * 50 * 10), dim = c(100, 50, 10))}\cr
#'  \code{names(dim(temperature)) <- c('longitude', 'latitude', 'time')}\cr
#'  \code{ArrayToNetCDF(list(temperature = temperature), file_path = 'example.nc')}
#'  }
#'  \item{
#'  Via the attribute 'variables' of each provided array: The arrays can 
#'  be provided with metadata in an attribute named 'variables', which is 
#'  expected to be a named list of named lists, where the names of the 
#'  container list are the names of the variables present in the provided 
#'  array, and where each sub-list contains metadata for each of the variables.
#'   The attribute names and values supported in the sub-lists must follow the 
#'  same format the package \code{ncdf4} uses to represent the NetCDF 
#'  file headers.\cr
#'  E.g:\cr
#'  \code{a <- array(1:400, dim = c(5, 10, 4, 2))}\cr
#'  \code{metadata <- list(tos = list(addOffset = 100,}\cr
#'  \code{                 scaleFact = 10,}\cr
#'  \code{                 dim = list(list(name = 'time',}\cr
#'  \code{                                 unlim = FALSE))))}\cr
#'  \code{attr(a, 'variables') <- metadata}\cr
#'  \code{names(dim(a)) <- c('lat', 'lon', 'time', 'var')}\cr
#'  \code{ArrayToNetCDF(a, 'tmp.nc')}
#'  }
#'}
#'The special dimension names are 'var'/'variable' and 'time'.\cr
#'If a dimension is named 'var' or 'variable', \code{ArrayToNetCDF} will 
#'interpret each array entry along such dimension corresponds to a separate 
#'new variable, hence will create a new variable inside the NetCDF file and 
#'will use it to store all the data in the provided array for the 
#'corresponding entry along the 'var'/'variable' dimension.\cr
#'If a dimension is named 'time', by default it will be interpreted and built 
#'as an unlimited dimension. The 'time' dimension must be the last dimension 
#'of the array (the right-most). If a 'var'/'variable' dimension is present, 
#'the 'time' dimension can be also placed on its left (i.e. the one before the 
#'last dimension). The default behaviour of creating the 'time' as unlimited 
#'dimension can be disabled by setting manually the attribute 
#'\code{unlim = FALSE}, as shown in the previous example.
#'
#'@param arrays One or a list of multidimensional data arrays. The list can be 
#'  provided with names, which will be interpreted as variable names. The 
#'  arrays can be provided with dimension names. The arrays can be provided 
#'  with metadata in the attribute 'variables' (read section Description for 
#'  details).
#'@param file_path Path and name of the NetCDF file to be created.
#'
#'@return This function returns NULL. 
#'@keywords datagen
#'@author History:\cr
#'  0.0 - 2017-01 (N. Manubens, \email{nicolau.manubens@bsc.es}) - Original code.
#'
#'@examples
#'  \dontrun{
#'# Minimal use case
#'ArrayToNetCDF(array(1:9, c(3, 3)), 'tmp.nc')
#'
#'# Works with arrays of any number of dimensions
#'ArrayToNetCDF(array(1:27, c(3, 3, 3)), 'tmp.nc')
#'
#'# Arrays can also be provided in [named] lists
#'ArrayToNetCDF(list(tos = array(1:27, c(3, 3, 3))), 'tmp.nc')
#'
#'# Or with dimension names
#'# 'var' dimension name will generate multiple variables in the 
#'# resulting NetCDF file
#'a <- array(1:27, dim = c(3, 3, 3))
#'names(dim(a)) <- c('lon', 'lat', 'var')
#'ArrayToNetCDF(a, 'tmp.nc')
#'
#'# 'variable' as dimension name will do the same
#'a <- array(1:27, dim = c(3, 3, 3))
#'names(dim(a)) <- c('lon', 'lat', 'variable')
#'ArrayToNetCDF(a, 'tmp.nc')
#'
#'# The 'time' dimension will be built as unlimited dimension, by default
#'a <- array(1:1600, dim = c(10, 20, 4, 2))
#'names(dim(a)) <- c('lat', 'lon', 'time', 'var')
#'ArrayToNetCDF(a, 'tmp.nc')
#'
#'# Putting the 'time' dimension in a position which is not the last, or the one
#'# right before 'var'/'variable' will crash. Unlimited dimension must be in the
#'# last position
#'a <- array(1:1600, dim = c(10, 20, 4, 2))
#'names(dim(a)) <- c('time', 'lat', 'lon', 'var')
#'ArrayToNetCDF(a, 'tmp.nc')
#'a <- array(1:1600, dim = c(10, 20, 4, 2))
#'names(dim(a)) <- c('lat', 'time', 'lon', 'var')
#'ArrayToNetCDF(a, 'tmp.nc')
#'
#'# The dimension 'var'/'variable' can be in any position and can have any length
#'a <- array(1:1600, dim = c(10, 20, 4, 2))
#'names(dim(a)) <- c('lat', 'var', 'lon', 'time')
#'ArrayToNetCDF(a, 'tmp.nc')
#'
#'# Multiple arrays can be provided in a list
#'a <- array(1:400, dim = c(5, 10, 4, 2))
#'names(dim(a)) <- c('lat', 'lon', 'time', 'var')
#'ArrayToNetCDF(list(a, a), 'tmp.nc')
#'
#'# If no dimension names are given to an array, new names will be automatically
#'# generated
#'a <- array(1:400, dim = c(5, 10, 4, 2))
#'b <- array(1:400, dim = c(5, 11, 4, 2))
#'names(dim(a)) <- c('lat', 'lon', 'time', 'var')
#'ArrayToNetCDF(list(a, b), 'tmp.nc')
#'
#'# If two arrays use a same dimension but their lengths differ, the function 
#'# will crash
#'a <- array(1:400, dim = c(5, 10, 4, 2))
#'b <- array(1:400, dim = c(5, 11, 4, 2))
#'names(dim(a)) <- c('lat', 'lon', 'time', 'var')
#'names(dim(b)) <- c('lat', 'lon', 'time', 'var')
#'ArrayToNetCDF(list(a, b), 'tmp.nc')
#'
#'# Metadata can be provided for each variable in each array, via the
#'# attribute 'variables'. In this example the metadata is empty.
#'a <- array(1:400, dim = c(5, 10, 4, 2))
#'metadata <- list(
#'              tos = list(),
#'              tas = list()
#'            )
#'attr(a, 'variables') <- metadata
#'names(dim(a)) <- c('lat', 'lon', 'time', 'var')
#'ArrayToNetCDF(a, 'tmp.nc')
#'
#'# Variable names can be manually specified
#'a <- array(1:400, dim = c(5, 10, 4, 2))
#'metadata <- list(
#'              tos = list(name = 'name1'),
#'              tas = list(name = 'name2')
#'            )
#'attr(a, 'variables') <- metadata
#'names(dim(a)) <- c('lat', 'lon', 'time', 'var')
#'ArrayToNetCDF(a, 'tmp.nc')
#'
#'# Units can be specified
#'a <- array(1:400, dim = c(5, 10, 4, 2))
#'metadata <- list(
#'              tos = list(units = 'K'),
#'              tas = list(units = 'K')
#'            )
#'attr(a, 'variables') <- metadata
#'names(dim(a)) <- c('lat', 'lon', 'time', 'var')
#'ArrayToNetCDF(a, 'tmp.nc')
#'
#'# addOffset and scaleFactor can be specified
#'a <- array(1:400, dim = c(5, 10, 4, 2))
#'metadata <- list(
#'              tos = list(addOffset = 100,
#'                         scaleFact = 10),
#'              tas = list(addOffset = 100,
#'                         scaleFact = 10)
#'            )
#'attr(a, 'variables') <- metadata
#'names(dim(a)) <- c('lat', 'lon', 'time', 'var')
#'ArrayToNetCDF(a, 'tmp.nc')
#'
#'# Unlimited dimensions can be manually created
#'a <- array(1:400, dim = c(5, 10, 4, 2))
#'metadata <- list(
#'              tos = list(addOffset = 100,
#'                         scaleFact = 10,
#'                         dim = list(list(name = 'unlimited',
#'                                         unlim = TRUE))),
#'              tas = list(addOffset = 100,
#'                         scaleFact = 10,
#'                         dim = list(list(name = 'unlimited',
#'                                         unlim = TRUE)))
#'            )
#'attr(a, 'variables') <- metadata
#'names(dim(a)) <- c('lat', 'lon', 'unlimited', 'var')
#'ArrayToNetCDF(a, 'tmp.nc')
#'
#'# A 'time' dimension can be built without it necessarily being unlimited
#'a <- array(1:400, dim = c(5, 10, 4, 2))
#'metadata <- list(
#'              tos = list(addOffset = 100,
#'                         scaleFact = 10,
#'                         dim = list(list(name = 'time',
#'                                         unlim = FALSE))),
#'              tas = list(addOffset = 100,
#'                         scaleFact = 10,
#'                         dim = list(list(name = 'time',
#'                                         unlim = FALSE)))
#'            )
#'attr(a, 'variables') <- metadata
#'names(dim(a)) <- c('lat', 'lon', 'time', 'var')
#'ArrayToNetCDF(a, 'tmp.nc')
#'
#'# Multiple arrays with data for multiple variables can be saved into a 
#'# NetCDF file at once.
#'tos <- array(1:400, dim = c(5, 10, 4))
#'metadata <- list(tos = list(units = 'K'))
#'attr(tos, 'variables') <- metadata
#'names(dim(tos)) <- c('lat', 'lon', 'time')
#'lon <- seq(0, 360 - 360 / 10, length.out = 10)
#'dim(lon) <- length(lon)
#'metadata <- list(lon = list(units = 'degrees_east'))
#'attr(lon, 'variables') <- metadata
#'names(dim(lon)) <- 'lon'
#'lat <- seq(-90, 90, length.out = 5)
#'dim(lat) <- length(lat)
#'metadata <- list(lat = list(units = 'degrees_north'))
#'attr(lat, 'variables') <- metadata
#'names(dim(lat)) <- 'lat'
#'ArrayToNetCDF(list(tos, lon, lat), 'tmp.nc')
#'}
#'@import ncdf4
#'@export
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
        if (any(is.na(dim_names) | (sapply(dim_names, nchar) == 0))) {
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
            dim_info[['len']] <- as.integer(round(dim_info[['len']][1]))
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
            if (!(is.numeric(dim_info[['vals']]))) {
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
        coords <- strsplit(var_info[['coordinates']], ' ')[[1]]
        if (!(all(coords %in% sapply(defined_vars, '[[', 'name') | 
                  coords %in% sapply(defined_dims[which(sapply(defined_dims, '[[', 'create_dimvar'))], '[[', 'name')))) {
          coords <- coords[which(coords %in% sapply(defined_vars, '[[', 'name') |
                                 coords %in% sapply(defined_dims[which(sapply(defined_dims, '[[', 'create_dimvar'))], '[[', 'name'))]
          .warning("Some of the dimensions appearing in 'coordinates' have been removed because they point to undefined variables.")
        }
        ncatt_put(ncdf_object, defined_vars[[var_counter]]$name, 'coordinates', paste(coords, collapse = ' '))
      }
      attrs_to_skip <- which(names(var_info) %in% c('addOffset', 'scaleFact', 'coordinates', 'dim'))
      attrs_to_add <- names(var_info)
      if (length(attrs_to_skip) > 0) {
        attrs_to_add <- attrs_to_add[-attrs_to_skip]
      }
      for (attribute_name in attrs_to_add) {
        if (is.numeric(var_info[[attribute_name]]) ||
            is.character(var_info[[attribute_name]])) {
          ncatt_put(ncdf_object, defined_vars[[var_counter]]$name, attribute_name, var_info[[attribute_name]])
        }
      }
      var_counter <- var_counter + 1
    }
  }
  nc_close(ncdf_object)
  invisible(NULL)
}
