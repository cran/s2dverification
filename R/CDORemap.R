CDORemap <- function(data_array = NULL, lons, lats, grid, method, 
                     avoid_writes = TRUE, crop = TRUE,
                     force_remap = FALSE, write_dir = tempdir()) {  #, mask = NULL) {
  .isRegularVector <- function(x, tol = 0.1) {
    if (length(x) < 2) {
      #stop("The provided vector must be of length 2 or greater.")
      TRUE
    } else {
      spaces <- x[2:length(x)] - x[1:(length(x) - 1)]
      (sum(abs(spaces - mean(spaces)) > mean(spaces) / (1 / tol)) < 2)
    }
  }
  # Check parameters data_array, lons and lats.
  known_lon_names <- .KnownLonNames()
  known_lat_names <- .KnownLatNames()
  if (!is.numeric(lons) || !is.numeric(lats)) {
    stop("Expected numeric 'lons' and 'lats'.")
  }
  if (any(is.na(lons > 0))) {
    stop("Found invalid values in 'lons'.")
  }
  if (any(is.na(lats > 0))) {
    stop("Found invalid values in 'lats'.")
  }
  if (is.null(dim(lons))) {
    dim(lons) <- length(lons)
  }
  if (is.null(dim(lats))) {
    dim(lats) <- length(lats)
  }
  if (length(dim(lons)) > 2 || length(dim(lats)) > 2) {
    stop("'lons' and 'lats' can only have up to 2 dimensions.")
  }
  if (length(dim(lons)) != length(dim(lats))) {
    stop("'lons' and 'lats' must have the same number of dimensions.")
  }
  if (length(dim(lons)) == 2 && !all(dim(lons) == dim(lats))) {
    stop("'lons' and 'lats' must have the same dimension sizes.")
  }
  return_array <- TRUE
  if (is.null(data_array)) {
    return_array <- FALSE
    if (length(dim(lons)) == 1) {
      array_dims <- c(length(lats), length(lons))
      names(array_dims) <- c('lat', 'lon')
    } else {
      array_dims <- dim(lons)
      names(array_dims) <- c('j', 'i')
    }
    data_array <- array(NA, array_dims)
  }
  if (!(is.logical(data_array) || is.numeric(data_array)) || !is.array(data_array)) {
    stop("Parameter 'data_array' must be a numeric array.")
  }
  if (is.null(names(dim(data_array)))) {
    stop("Parameter 'data_array' must have named dimensions.")
  }
  lon_dim <- which(known_lon_names %in% names(dim(data_array)))
  if (length(lon_dim) < 1) {
    stop("Could not find a known longitude dimension name in the provided 'data_array'.")
  }
  if (length(lon_dim) > 1) {
    stop("Found more than one known longitude dimension names in the provided 'data_array'.")
  }
  lon_dim <- known_lon_names[lon_dim]
  lat_dim <- which(known_lat_names %in% names(dim(data_array)))
  if (length(lat_dim) < 1) {
    stop("Could not find a known latitude dimension name in the provided 'data_array'.")
  }
  if (length(lat_dim) > 1) {
    stop("Found more than one known latitude dimension name in the provided 'data_array'.")
  }
  lat_dim <- known_lat_names[lat_dim]
  if (is.null(names(dim(lons)))) {
    if (length(dim(lons)) == 1) {
      names(dim(lons)) <- lon_dim
    } else {
      stop("Parameter 'lons' must be provided with dimension names.")
    }
  } else {
    if (!(lon_dim %in% names(dim(lons)))) {
      stop("Parameter 'lon' must have the same longitude dimension name as the 'data_array'.")
    }
    if (length(dim(lons)) > 1 && !(lat_dim %in% names(dim(lons)))) {
      stop("Parameter 'lon' must have the same latitude dimension name as the 'data_array'.")
    }
  }
  if (is.null(names(dim(lats)))) {
    if (length(dim(lats)) > 1) {
      stop("Parameter 'lats' must be provided with dimension names.")
    }
  } else {
    if (!(lat_dim %in% names(dim(lats)))) {
      stop("Parameter 'lat' must have the same latitude dimension name as the 'data_array'.")
    }
    if (length(dim(lats)) > 1 && !(lon_dim %in% names(dim(lats)))) {
      stop("Parameter 'lat' must have the same longitude dimension name as the 'data_array'.")
    }
  }
  lons_attr_bk <- attributes(lons)
  if (is.null(lons_attr_bk)) {
    lons_attr_bk <- list()
  }
  lats_attr_bk <- attributes(lats)
  if (is.null(lats_attr_bk)) {
    lats_attr_bk <- list()
  }
  if (length(attr(lons, 'variables')) == 0) {
    new_metadata <- list(list())
    if (length(dim(lons)) == 1) {
      names(new_metadata) <- lon_dim
    } else {
      names(new_metadata) <- paste0(lon_dim, '_var')
    }
    attr(lons, 'variables') <- new_metadata
  }
  if (!('units' %in% names(attr(lons, 'variables')[[1]]))) {
    new_metadata <- attr(lons, 'variables')
    #names(new_metadata)[1] <- lon_dim
    new_metadata[[1]][['units']] <- 'degrees_east'
    attr(lons, 'variables') <- new_metadata
  }
  if (length(attr(lats, 'variables')) == 0) {
    new_metadata <- list(list())
    if (length(dim(lats)) == 1) {
      names(new_metadata) <- lat_dim
    } else {
      names(new_metadata) <- paste0(lat_dim, '_var')
    }
    attr(lats, 'variables') <- new_metadata
  }
  if (!('units' %in% names(attr(lats, 'variables')[[1]]))) {
    new_metadata <- attr(lats, 'variables')
    #names(new_metadata)[1] <- lat_dim
    new_metadata[[1]][['units']] <- 'degrees_north'
    attr(lats, 'variables') <- new_metadata
  }
  # Check grid.
  if (!is.character(grid)) {
    stop("Parameter 'grid' must be a character string specifying a ",
         "target CDO grid, 'rXxY' or 'tRESgrid', or a path to another ",
         "NetCDF file.")
  }
  if (grepl('^r[0-9]{1,}x[0-9]{1,}$', grid)) {
    grid_type <- 'regular'
    grid_lons <- as.numeric(strsplit(strsplit(grid, 'x')[[1]][1], 'r')[[1]][2])
    grid_lats <- as.numeric(strsplit(grid, 'x')[[1]][2])
  } else if (grepl('^t[0-9]{1,}grid$', grid)) {
    grid_type <- 'gaussian'
    grid_t <- as.numeric(strsplit(strsplit(grid, 'grid')[[1]][1], 't')[[1]][2])
    grid_size <- .t2nlatlon(grid_t)
    grid_lons <- grid_size[2]
    grid_lats <- grid_size[1]
  } else {
    grid_type <- 'custom'
  }
  # Check method.
  if (method %in% c('bil', 'bilinear')) {
    method <- 'bil'
  } else if (method %in% c('bic', 'bicubic')) {
    method <- 'bic'
  } else if (method %in% c('con', 'conservative')) {
    method <- 'con'
  } else if (method %in% c('dis', 'distance-weighted')) {
    method <- 'dis'
  } else {
    stop("Unsupported CDO remap method. 'bilinear', 'bicubic', 'conservative' or 'distance-weighted' supported only.")
  }
  # Check avoid_writes
  if (!is.logical(avoid_writes)) {
    stop("Parameter 'avoid_writes' must be a logical value.")
  }
  # Check crop
  crop_tight <- FALSE
  if (is.character(crop)) {
    if (crop == 'tight') {
      crop_tight <- TRUE
    } else if (crop != 'preserve') {
      stop("Parameter 'crop' can only take the values 'tight' or 'preserve' if specified as a character string.")
    }
    crop <- TRUE
  }
  if (is.logical(crop)) {
    if (crop) {
      if (length(lons) == 1 || length(lats) == 1) {
        stop("CDORemap cannot remap if crop = TRUE and values for only one ",
             "longitude or one latitude are provided. Either a) provide ",
             "values for more than one longitude/latitude, b) explicitly ",
             "specify the crop limits in the parameter crop, or c) set ",
             "crop = FALSE.")
      }
      if (crop_tight) {
        lon_extremes <- c(min(lons), max(lons))
        lat_extremes <- c(min(lats), max(lats))
      } else {
        # Here we are trying to look for the extreme lons and lats in the data.
        # Not the centers of the extreme cells, but the borders of the extreme cells.
###---
        if (length(dim(lons)) == 1) {
          tmp_lon <- lons
        } else {
          min_pos <- which(lons == min(lons), arr.ind = TRUE)[1, ]
          tmp_lon <- Subset(lons, lat_dim, min_pos[which(names(dim(lons)) == lat_dim)], drop = 'selected')
        }
        i <- 1:length(tmp_lon)
        lon_model <- lm(tmp_lon ~ poly(i, 3))
        lon_extremes <- c(NA, NA)
        left_is_min <- FALSE
        right_is_max <- FALSE
        if (which.min(tmp_lon) == 1) {
          left_is_min <- TRUE
          prev_lon <- predict(lon_model, data.frame(i = 0))
          first_lon_cell_width <- (tmp_lon[1] - prev_lon)
          lon_extremes[1] <- tmp_lon[1] - first_lon_cell_width / 2
        } else {
          lon_extremes[1] <- min(tmp_lon)
        }
        if (which.max(tmp_lon) == length(tmp_lon)) {
          right_is_max <- TRUE
          next_lon <- predict(lon_model, data.frame(i = length(tmp_lon) + 1))
          last_lon_cell_width <- (next_lon - tmp_lon[length(tmp_lon)])
          lon_extremes[2] <- tmp_lon[length(tmp_lon)] + last_lon_cell_width / 2
        } else {
          lon_extremes[2] <- max(tmp_lon)
        }
        # Adjust the crop window if possible in order to keep lons from 0 to 360 
        # or from -180 to 180 when the extremes of the cropped window are contiguous.
        if (right_is_max) {
          if (lon_extremes[1] < -180) {
            if (!((lon_extremes[2] < 180) && !((180 - lon_extremes[2]) <= last_lon_cell_width / 2))) {
              lon_extremes[1] <- -180
              lon_extremes[2] <- 180
            }
          } else if (lon_extremes[1] < 0) {
            if (!((lon_extremes[2] < 360) && !((360 - lon_extremes[2]) <= last_lon_cell_width / 2))) {
              lon_extremes[1] <- 0
              lon_extremes[2] <- 360
            }
          }
        } 
        if (left_is_min) {
          if (lon_extremes[2] > 360) {
            if (!((lon_extremes[1] > 0) && !(lon_extremes[1] <= first_lon_cell_width / 2))) {
              lon_extremes[1] <- 0
              lon_extremes[2] <- 360
            }
          } else if (lon_extremes[2] > 180) {
            if (!((lon_extremes[1] > -180) && !((180 + lon_extremes[1]) <= first_lon_cell_width / 2))) {
              lon_extremes[1] <- -180
              lon_extremes[2] <- 180
            }
          }
        } 
##      lon_extremes <- signif(lon_extremes, 5)
##      lon_extremes <- lon_extremes + 0.00001
###---
        if (length(dim(lats)) == 1) {
          tmp_lat <- lats
        } else {
          min_pos <- which(lats == min(lats), arr.ind = TRUE)[1, ]
          tmp_lat <- Subset(lats, lon_dim, min_pos[which(names(dim(lats)) == lon_dim)], drop = 'selected')
        }
        i <- 1:length(tmp_lat)
        lat_model <- lm(tmp_lat ~ poly(i, 3))
        lat_extremes <- c(NA, NA)
        if (which.min(tmp_lat) == 1) {
          prev_lat <- predict(lat_model, data.frame(i = 0))
          lat_extremes[1] <- tmp_lat[1] - (tmp_lat[1] - prev_lat) / 2
        } else {
          lat_extremes[1] <- min(tmp_lat)
        }
        if (which.max(tmp_lat) == length(tmp_lat)) {
          next_lat <- predict(lat_model, data.frame(i = length(tmp_lat) + 1))
          lat_extremes[2] <- tmp_lat[length(tmp_lat)] + (next_lat - tmp_lat[length(tmp_lat)]) / 2
        } else {
          lat_extremes[2] <- max(tmp_lat)
        }
##      lat_extremes <- signif(lat_extremes, 5)
        # Adjust crop window
        if (lat_extremes[1] < -90) {
          lat_extremes[1] <- -90
        } else if (lat_extremes[1] > 90) {
          lat_extremes[1] <- 90
        }
        if (lat_extremes[2] < -90) {
          lat_extremes[2] <- -90
        } else if (lat_extremes[2] > 90) {
          lat_extremes[2] <- 90
        }
###---
      }
    }
  } else if (is.numeric(crop)) {
    if (length(crop) != 4) {
      stop("Paramrter 'crop' must be a logical value or a numeric vector of length 4: c(western border, eastern border, southern border, northern border.")
    } else {
      lon_extremes <- crop[1:2]
      lat_extremes <- crop[3:4]
      crop <- TRUE
    }
  } else {
    stop("Parameter 'crop' must be a logical value or a numeric vector.")
  }
  # Check force_remap
  if (!is.logical(force_remap)) {
    stop("Parameter 'force_remap' must be a logical value.")
  }
  # Check write_dir
  if (!is.character(write_dir)) {
    stop("Parameter 'write_dir' must be a character string.")
  }
  if (!dir.exists(write_dir)) {
    stop("Parameter 'write_dir' must point to an existing directory.")
  }
#  if (!is.null(mask)) {
#    if (!is.numeric(mask) || !is.array(mask)) {
#      stop("Parameter 'mask' must be a numeric array.")
#    }
#    if (length(dim(mask)) != 2) {
#      stop("Parameter 'mask' must have two dimensions.")
#    }
#    if (is.null(names(dim(mask)))) {
#      if (dim(data_array)[lat_dim] == dim(data_array)[lon_dim]) {
#        stop("Cannot disambiguate which is the longitude dimension of ",
#             "the provided 'mask'. Provide it with dimension names.")
#      }
#      names(dim(mask)) <- c('', '')
#      found_lon_dim <- which(dim(mask) == dim(data_array)[lon_dim])
#      if (length(found_lon_dim) < 0) {
#        stop("The dimension sizes of the provided 'mask' do not match ",
#             "the spatial dimension sizes of the array to interpolate.")
#      } else {
#        names(dim(mask)[found_lon_dim]) <- lon_dim
#      }
#      found_lat_dim <- which(dim(mask) == dim(data_array)[lat_dim])
#      if (length(found_lat_dim) < 0) {
#        stop("The dimension sizes of the provided 'mask' do not match ",
#             "the spatial dimension sizes of the array to interpolate.")
#      } else {
#        names(dim(mask)[found_lat_dim]) <- lat_dim
#      }
#    }
#    lon_position <- which(names(dim(data_array)) == lon_dim)
#    lat_position <- which(names(dim(data_array)) == lat_dim)
#    if (lon_position > lat_position) {
#      if (names(dim(mask))[1] == lon_dim) {
#        mask <- t(mask)
#      }
#    } else {
#      if (names(dim(mask))[1] == lat_dim) {
#        mask <- t(mask)
#      }
#    }
#    ## TODO: Apply mask!!! Preserve attributes
#  }
  # Check if interpolation can be skipped.
  interpolation_needed <- TRUE
  if (!force_remap) {
    if (!(grid_type == 'custom')) {
      if (length(lons) == grid_lons && length(lats) == grid_lats) {
        if (grid_type == 'regular') {
          if (.isRegularVector(lons) && .isRegularVector(lats)) {
            interpolation_needed <- FALSE
          }
        } else if (grid_type == 'gaussian') {
          # TODO: improve this check. Gaussian quadrature should be used.
          if (.isRegularVector(lons) && !.isRegularVector(lats)) {
            interpolation_needed <- FALSE
          }
        }
      }
    }
  }
  found_lons <- lons
  found_lats <- lats
  if (interpolation_needed) {
    if (nchar(Sys.which('cdo')[1]) < 1) {
      stop("CDO must be installed in order to use the .CDORemap.")
    }
    # CDO takes arrays of 3 dimensions or 4 if one of them is unlimited.
    # The unlimited dimension can only be the left-most (right-most in R).
    # There are no restrictions for the dimension names or variable names.
    # The longitude and latitude are detected by their units.
    # There are no restrictions for the order of the limited dimensions.
    # The longitude/latitude variables and dimensions must have the same name.
    # The procedure consists in:
    # - take out the array metadata
    # - be aware of var dimension (replacing the dimension names would do).
    # - take arrays of 4 dimensions always if possible
    # - make the last dimension unlimited when saving to netcdf
    # - if the last dimension is lon or lat, either reorder the array and 
    #   then reorder back or iterate over the dimensions at the right
    #   side of lon AND lat.
    # If the input array has more than 4 dimensions, it is needed to
    # run CDO on each sub-array of 4 dimensions because it can handle
    # only up to 4 dimensions. The shortest dimensions are chosen to 
    # iterate over.
    is_irregular <- FALSE
    if (length(dim(lats)) > 1 && length(dim(lons)) > 1) {
      is_irregular <- TRUE
    }
    attribute_backup <- attributes(data_array)
    other_dims <- which(!(names(dim(data_array)) %in% c(lon_dim, lat_dim)))
    permutation <- NULL
    unlimited_dim <- NULL
    dims_to_iterate <- NULL
    total_slices <- 1
    other_dims_per_chunk <- ifelse(is_irregular, 1, 2)  # 4 (the maximum accepted by CDO) - 2 (lon, lat) = 2.
    if (length(other_dims) > 1 || (length(other_dims) > 0 && (is_irregular))) {
      if (!(length(dim(data_array)) %in% other_dims)) {
        if (avoid_writes || is_irregular) {
          dim_to_move <- max(other_dims)
          permutation <- (1:length(dim(data_array)))[-dim_to_move]
          permutation <- c(permutation, dim_to_move)
          permutation_back <- 1:length(dim(data_array))
          permutation_back[dim_to_move] <- length(dim(data_array))
          permutation_back[length(dim(data_array))] <- dim_to_move
          dim_backup <- dim(data_array)
          data_array <- aperm(data_array, permutation)
          dim(data_array) <- dim_backup[permutation]
          other_dims <- which(!(names(dim(data_array)) %in% c(lon_dim, lat_dim)))
        } else {
          # We allow only lon, lat and 1 more dimension per chunk, so 
          # CDO has no restrictions in the order.
          other_dims_per_chunk <- 1
        }
      }
      other_dims_ordered_by_size <- other_dims[sort(dim(data_array)[other_dims], index.return = TRUE)$ix]
      dims_to_iterate <- sort(head(other_dims_ordered_by_size, length(other_dims) - other_dims_per_chunk))
      if (length(dims_to_iterate) == 0) {
        dims_to_iterate <- NULL
      } else {
        slices_to_iterate <- array(1:prod(dim(data_array)[dims_to_iterate]), 
                                    dim(data_array)[dims_to_iterate])
        total_slices <- prod(dim(slices_to_iterate))
      }
      if ((other_dims_per_chunk > 1) || (other_dims_per_chunk > 0 && is_irregular)) {
        unlimited_dim <- tail(sort(tail(other_dims_ordered_by_size, other_dims_per_chunk)), 1)
      }
    }

    result_array <- NULL
    lon_pos <- which(names(dim(data_array)) == lon_dim)
    lat_pos <- which(names(dim(data_array)) == lat_dim)
    dim_backup <- dim(data_array)
    attributes(data_array) <- NULL
    dim(data_array) <- dim_backup
    names(dim(data_array)) <- paste0('dim', 1:length(dim(data_array)))
    names(dim(data_array))[c(lon_pos, lat_pos)] <- c(lon_dim, lat_dim)
    if (!is.null(unlimited_dim)) {
      # This will make ArrayToNetCDF create this dim as unlimited.
      names(dim(data_array))[unlimited_dim] <- 'time'
    }
    if (length(dim(lons)) == 1) {
      names(dim(lons)) <- lon_dim
    }
    if (length(dim(lats)) == 1) {
      names(dim(lats)) <- lat_dim
    }
    if (length(dim(lons)) > 1) {
      lon_var_name <- paste0(lon_dim, '_var')
    } else {
      lon_var_name <- lon_dim
    }
    if (length(dim(lats)) > 1) {
      lat_var_name <- paste0(lat_dim, '_var')
    } else {
      lat_var_name <- lat_dim
    }
    if (is_irregular) {
      metadata <- list(list(coordinates = paste(lon_var_name, lat_var_name)))
      names(metadata) <- 'var'
      attr(data_array, 'variables') <- metadata
    }
    names(attr(lons, 'variables')) <- lon_var_name
    names(attr(lats, 'variables')) <- lat_var_name
    for (i in 1:total_slices) {
      tmp_file <- tempfile('R_CDORemap_', write_dir, fileext = '.nc')
      tmp_file2 <- tempfile('R_CDORemap_', write_dir, fileext = '.nc')
      if (!is.null(dims_to_iterate)) {
        slice_indices <- which(slices_to_iterate == i, arr.ind = TRUE)
        subset <- Subset(data_array, dims_to_iterate, as.list(slice_indices), drop = 'selected')
        # Make sure subset goes along with metadata
        ArrayToNetCDF(setNames(list(subset, lons, lats), c('var', lon_var_name, lat_var_name)), tmp_file)
      } else {
        ArrayToNetCDF(setNames(list(data_array, lons, lats), c('var', lon_var_name, lat_var_name)), tmp_file)
      }
      sellonlatbox <- ''
      if (crop) {
        sellonlatbox <- paste0('sellonlatbox,', lon_extremes[1], ',', lon_extremes[2], 
                                           ',', lat_extremes[1], ',', lat_extremes[2], ' -')
      }
      err <- try({
## TODO: Here add sellonlatbox. Also check constantin's issue, may contain hint. Also search if possible to crop without
        system(paste0("cdo -s ", sellonlatbox, "remap", method, ",", grid, " ", tmp_file, " ", tmp_file2))
      })
      file.remove(tmp_file)
      if (('try-error' %in% class(err)) || err > 0) {
        stop("CDO remap failed.")
      }
      ncdf_remapped <- nc_open(tmp_file2)
      found_dim_names <- sapply(ncdf_remapped$var$var$dim, '[[', 'name')
      found_lon_dim <- found_dim_names[which(found_dim_names %in% .KnownLonNames())[1]]
      found_lat_dim <- found_dim_names[which(found_dim_names %in% .KnownLatNames())[1]]
      found_lon_dim_size <- length(ncdf_remapped$dim[[found_lon_dim]]$vals)
      found_lat_dim_size <- length(ncdf_remapped$dim[[found_lat_dim]]$vals)
      found_lons <- ncvar_get(ncdf_remapped, 'lon', collapse_degen = FALSE)
      found_lats <- ncvar_get(ncdf_remapped, 'lat', collapse_degen = FALSE)
      if (length(dim(found_lons)) > 1) {
        if (found_lon_dim < found_lat_dim) {
          names(dim(found_lons)) <- c(found_lon_dim, found_lat_dim)
        } else {
          names(dim(found_lons)) <- c(found_lat_dim, found_lon_dim)
        }
      } else {
        names(dim(found_lons)) <- found_lon_dim
      }
      if (length(dim(found_lats)) > 1) {
        if (found_lon_dim < found_lat_dim) {
          names(dim(found_lats)) <- c(found_lon_dim, found_lat_dim)
        } else {
          names(dim(found_lats)) <- c(found_lat_dim, found_lon_dim)
        }
      } else {
        names(dim(found_lats)) <- found_lat_dim
      }
      if (!is.null(dims_to_iterate)) {
        if (is.null(result_array)) {
          if (return_array) {
            new_dims <- dim(data_array)
            new_dims[c(lon_dim, lat_dim)] <- c(found_lon_dim_size, found_lat_dim_size)
            result_array <- array(dim = new_dims)
            store_indices <- as.list(rep(TRUE, length(dim(result_array))))
          }
        }
        if (return_array) {
          store_indices[dims_to_iterate] <- as.list(slice_indices)
          result_array <- do.call('[<-', c(list(x = result_array), store_indices, 
                                           list(value = ncvar_get(ncdf_remapped, 'var', collapse_degen = FALSE))))
        }
      } else {
        new_dims <- dim(data_array)
        new_dims[c(lon_dim, lat_dim)] <- c(found_lon_dim_size, found_lat_dim_size)
        result_array <- ncvar_get(ncdf_remapped, 'var', collapse_degen = FALSE)
        names(dim(result_array)) <- names(new_dims)
      }
      nc_close(ncdf_remapped)
      file.remove(tmp_file2)
    }
    if (!is.null(permutation)) {
      dim_backup <- dim(result_array)
      result_array <- aperm(result_array, permutation_back)
      dim(result_array) <- dim_backup[permutation_back]
    }
    # Now restore the metadata
    result_is_irregular <- FALSE
    if (length(dim(found_lats)) > 1 && length(dim(found_lons)) > 1) {
      result_is_irregular <- TRUE
    }
    attribute_backup[['dim']][which(names(dim(result_array)) == lon_dim)] <- dim(result_array)[lon_dim]
    attribute_backup[['dim']][which(names(dim(result_array)) == lat_dim)] <- dim(result_array)[lat_dim]
    new_lon_name <- names(dim(found_lons))[which(names(dim(found_lons)) %in% .KnownLonNames())]
    new_lat_name <- names(dim(found_lats))[which(names(dim(found_lats)) %in% .KnownLatNames())]
    names(attribute_backup[['dim']])[which(names(dim(result_array)) == lon_dim)] <- new_lon_name
    names(attribute_backup[['dim']])[which(names(dim(result_array)) == lat_dim)] <- new_lat_name
    if (!is.null(attribute_backup[['variables']]) && (length(attribute_backup[['variables']]) > 0)) {
      for (var in 1:length(attribute_backup[['variables']])) {
        if (length(attribute_backup[['variables']][[var]][['dim']]) > 0) {
          for (dim in 1:length(attribute_backup[['variables']][[var]][['dim']])) {
            dim_name <- NULL
            if ('name' %in% names(attribute_backup[['variables']][[var]][['dim']][[dim]])) {
              dim_name <- attribute_backup[['variables']][[var]][['dim']][[dim]][['name']]
              if (dim_name %in% c(lon_dim, lat_dim)) {
                if (dim_name == lon_dim) {
                  attribute_backup[['variables']][[var]][['dim']][[dim]][['name']] <- new_lon_name
                } else {
                  attribute_backup[['variables']][[var]][['dim']][[dim]][['name']] <- new_lat_name
                }
              }
            } else if (!is.null(names(attribute_backup[['variables']][[var]][['dim']]))) {
              dim_name <- names(attribute_backup[['variables']][[var]][['dim']])[dim]
              if (dim_name %in% c(lon_dim, lat_dim)) {
                if (dim_name == lon_dim) {
                  names(attribute_backup[['variables']][[var]][['dim']])[which(names(attribute_backup[['variables']][[var]][['dim']]) == lon_dim)] <- new_lon_name
                } else {
                  names(attribute_backup[['variables']][[var]][['dim']])[which(names(attribute_backup[['variables']][[var]][['dim']]) == lat_dim)] <- new_lat_name
                }
              }
            }
            if (!is.null(dim_name)) {
              if (dim_name %in% c(lon_dim, lat_dim)) {
                if (dim_name == lon_dim) {
                  new_vals <- found_lons[TRUE]
                } else if (dim_name == lat_dim) {
                  new_vals <- found_lats[TRUE]
                }
                if (!is.null(attribute_backup[['variables']][[var]][['dim']][[dim]][['len']])) {
                  attribute_backup[['variables']][[var]][['dim']][[dim]][['len']] <- length(new_vals)
                }
                if (!is.null(attribute_backup[['variables']][[var]][['dim']][[dim]][['vals']])) {
                  if (!result_is_irregular) {
                    attribute_backup[['variables']][[var]][['dim']][[dim]][['vals']] <- new_vals
                  } else {
                    attribute_backup[['variables']][[var]][['dim']][[dim]][['vals']] <- 1:length(new_vals)
                  }
                }
              }
            }
          }
        }
        if (!is_irregular && result_is_irregular) {
          attribute_backup[['coordinates']] <- paste(lon_var_name, lat_var_name)
        } else if (is_irregular && !result_is_irregular) {
          attribute_backup[['coordinates']] <- NULL
        }
      }
    }
    attributes(result_array) <- attribute_backup
    lons_attr_bk[['dim']] <- dim(found_lons)
    if (!is.null(lons_attr_bk[['variables']]) && (length(lons_attr_bk[['variables']]) > 0)) {
      for (var in 1:length(lons_attr_bk[['variables']])) {
        if (length(lons_attr_bk[['variables']][[var]][['dim']]) > 0) {
          dims_to_remove <- NULL
          for (dim in 1:length(lons_attr_bk[['variables']][[var]][['dim']])) {
            dim_name <- NULL
            if ('name' %in% names(lons_attr_bk[['variables']][[var]][['dim']][[dim]])) {
              dim_name <- lons_attr_bk[['variables']][[var]][['dim']][[dim]][['name']]
              if (dim_name %in% c(lon_dim, lat_dim)) {
                if (dim_name == lon_dim) {
                  lons_attr_bk[['variables']][[var]][['dim']][[dim]][['name']] <- new_lon_name
                } else {
                  lons_attr_bk[['variables']][[var]][['dim']][[dim]][['name']] <- new_lat_name
                }
              }
            } else if (!is.null(names(lons_attr_bk[['variables']][[var]][['dim']]))) {
              dim_name <- names(lons_attr_bk[['variables']][[var]][['dim']])[dim]
              if (dim_name %in% c(lon_dim, lat_dim)) {
                if (dim_name == lon_dim) {
                  names(lons_attr_bk[['variables']][[var]][['dim']])[which(names(lons_attr_bk[['variables']][[var]][['dim']]) == lon_dim)] <- new_lon_name
                } else {
                  names(lons_attr_bk[['variables']][[var]][['dim']])[which(names(lons_attr_bk[['variables']][[var]][['dim']]) == lat_dim)] <- new_lat_name
                }
              }
            }
            if (!is.null(dim_name)) {
              if (dim_name %in% c(lon_dim, lat_dim)) {
                if (dim_name == lon_dim) {
                  new_vals <- found_lons[TRUE]
                } else if (dim_name == lat_dim) {
                  new_vals <- found_lats[TRUE]
                  if (!result_is_irregular) {
                    dims_to_remove <- c(dims_to_remove, dim)
                  }
                }
                if (!is.null(lons_attr_bk[['variables']][[var]][['dim']][[dim]][['len']])) {
                  lons_attr_bk[['variables']][[var]][['dim']][[dim]][['len']] <- length(new_vals)
                }
                if (!is.null(lons_attr_bk[['variables']][[var]][['dim']][[dim]][['vals']])) {
                  if (!result_is_irregular) {
                    lons_attr_bk[['variables']][[var]][['dim']][[dim]][['vals']] <- new_vals
                  } else {
                    lons_attr_bk[['variables']][[var]][['dim']][[dim]][['vals']] <- 1:length(new_vals)
                  }
                }
              }
            }
          }
          if (length(dims_to_remove) > 1) {
            lons_attr_bk[['variables']][[var]][['dim']] <- lons_attr_bk[['variables']][[var]][['dim']][[-dims_to_remove]]
          }
        }
      }
      names(lons_attr_bk[['variables']])[1] <- lon_var_name
      lons_attr_bk[['variables']][[1]][['units']] <- 'degrees_east'
    }
    attributes(found_lons) <- lons_attr_bk
    lats_attr_bk[['dim']] <- dim(found_lats)
    if (!is.null(lats_attr_bk[['variables']]) && (length(lats_attr_bk[['variables']]) > 0)) {
      for (var in 1:length(lats_attr_bk[['variables']])) {
        if (length(lats_attr_bk[['variables']][[var]][['dim']]) > 0) {
          dims_to_remove <- NULL
          for (dim in 1:length(lats_attr_bk[['variables']][[var]][['dim']])) {
            dim_name <- NULL
            if ('name' %in% names(lats_attr_bk[['variables']][[var]][['dim']][[dim]])) {
              dim_name <- lats_attr_bk[['variables']][[var]][['dim']][[dim]][['name']]
              if (dim_name %in% c(lon_dim, lat_dim)) {
                if (dim_name == lon_dim) {
                  lons_attr_bk[['variables']][[var]][['dim']][[dim]][['name']] <- new_lon_name
                } else {
                  lons_attr_bk[['variables']][[var]][['dim']][[dim]][['name']] <- new_lat_name
                }
              }
            } else if (!is.null(names(lats_attr_bk[['variables']][[var]][['dim']]))) {
              dim_name <- names(lats_attr_bk[['variables']][[var]][['dim']])[dim]
              if (dim_name %in% c(lon_dim, lat_dim)) {
                if (dim_name == lon_dim) {
                  names(lats_attr_bk[['variables']][[var]][['dim']])[which(names(lats_attr_bk[['variables']][[var]][['dim']]) == lon_dim)] <- new_lon_name
                } else {
                  names(lats_attr_bk[['variables']][[var]][['dim']])[which(names(lats_attr_bk[['variables']][[var]][['dim']]) == lat_dim)] <- new_lat_name
                }
              }
            }
            if (!is.null(dim_name)) {
              if (dim_name %in% c(lon_dim, lat_dim)) {
                if (dim_name == lon_dim) {
                  new_vals <- found_lons[TRUE]
                  if (!result_is_irregular) {
                    dims_to_remove <- c(dims_to_remove, dim)
                  }
                } else if (dim_name == lat_dim) {
                  new_vals <- found_lats[TRUE]
                }
                if (!is.null(lats_attr_bk[['variables']][[var]][['dim']][[dim]][['len']])) {
                  lats_attr_bk[['variables']][[var]][['dim']][[dim]][['len']] <- length(new_vals)
                }
                if (!is.null(lats_attr_bk[['variables']][[var]][['dim']][[dim]][['vals']])) {
                  if (!result_is_irregular) {
                    lats_attr_bk[['variables']][[var]][['dim']][[dim]][['vals']] <- new_vals
                  } else {
                    lats_attr_bk[['variables']][[var]][['dim']][[dim]][['vals']] <- 1:length(new_vals)
                  }
                }
              }
            }
          }
          if (length(dims_to_remove) > 1) {
            lats_attr_bk[['variables']][[var]][['dim']] <- lats_attr_bk[['variables']][[var]][['dim']][[-dims_to_remove]]
          }
        }
      }
      names(lats_attr_bk[['variables']])[1] <- lat_var_name
      lats_attr_bk[['variables']][[1]][['units']] <- 'degrees_north'
    }
    attributes(found_lats) <- lats_attr_bk
  }
  list(data_array = if (return_array) {
                      if (interpolation_needed) {
                        result_array
                      } else {
                        data_array
                      }
                    } else {
                      NULL
                    },
       lons = found_lons, lats = found_lats)
}
