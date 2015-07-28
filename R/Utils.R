## Function to tell if a regexpr() match is a complete match to a specified name
.IsFullMatch <- function(x, name) {
  ifelse(x > 0 && attributes(x)$match.length == nchar(name), TRUE, FALSE)
}

.ConfigReplaceVariablesInString <- function(string, replace_values, allow_undefined_key_vars = FALSE) {
  # This function replaces all the occurrences of a variable in a string by 
  # their corresponding string stored in the replace_values.
  if (length(strsplit(string, "\\$")[[1]]) > 1) {
    parts <- strsplit(string, "\\$")[[1]]
    output <- ""
    i <- 0
    for (part in parts) {
      if (i %% 2 == 0) {
        output <- paste(output, part, sep = "")
      } else {
        if (part %in% names(replace_values)) {
          output <- paste(output, .ConfigReplaceVariablesInString(replace_values[[part]], replace_values, allow_undefined_key_vars), sep = "")
        } else if (allow_undefined_key_vars && (part %in% c("STORE_FREQ", "START_DATE", "MEMBER_NUMBER", "YEAR", "MONTH"))) {
          output <- paste0(output, "$", part, "$")
        } else {
          stop(paste('Error: The variable $', part, '$ was not defined in the configuration file.', sep = ''))
        }
      }
      i <- i + 1
    }
    output
  } else {
    string
  }
}

.LoadDataFile <- function(work_piece, explore_dims = FALSE, silent = FALSE) {
  #suppressPackageStartupMessages({library(ncdf4)})
  #suppressPackageStartupMessages({library(bigmemory)})

  arrayIndex2VectorIndex <- function(indices, dims) {
    if (length(indices) > length(dims)) {
      stop("Error: indices do not match dimensions in arrayIndex2VectorIndex.")
    }
    position <- 1
    dims <- rev(dims)
    indices <- rev(indices)
    for (i in 1:length(indices)) {
      position <- position + (indices[i] - 1) * prod(dims[-c(1:i)])
    }
    position
  }

  file_found <- FALSE
  dims <- NULL

  filename <- work_piece[['filename']]
  namevar <- work_piece[['namevar']]
  output <- work_piece[['output']]
  # The names of all data files in the directory of the repository that match 
  # the pattern are obtained.
  files <- system(paste('ls', filename), intern = TRUE, ignore.stderr = TRUE)
  
  # If we don't find any, we set the flag 'file_found' with a false value.
  if (length(files) > 0) {
    # The first file that matches the pattern is chosen and read.
    filename <- files[length(files)]
    file_found <- TRUE

    if (!silent) {
      if (explore_dims) {
        cat(paste("* Exploring dimensions...", filename, '\n'))
      }
      ##} else {
      ##  cat(paste("* Reading & processing data...", filename, '\n'))
      ##}
    }

    # The data is regridded if it corresponds to an atmospheric variable. When
    # the chosen output type is 'areave' the data is not regridded to not 
    # waste computing time unless the user specified a common grid.
    if (work_piece[['is_2d_var']] && !is.null(work_piece[['grid']])) {
      ## TODO: Here we should check if the file is already in the
      ## target grid and avoid useless interpolations.
      filein <- tempfile(pattern = "loadRegridded", fileext = ".nc")
      system(paste("cdo ", work_piece[['remap']], ",", work_piece[['grid']], " -selname,", namevar, " ", 
                   filename, " ", filein, " 2>/dev/null", sep = ""))
      system(paste0("ncrename -v lat,", work_piece[['dimnames']][['latitudes']], 
                    " -d lat,", work_piece[['dimnames']][['latitudes']], 
                    " -v lon,", work_piece[['dimnames']][['longitudes']],
                    " -d lon,", work_piece[['dimnames']][['longitudes']],
                    " -v .lev,", work_piece[['dimnames']][['members']],
                    " -d .lev,", work_piece[['dimnames']][['members']],
                    " ", filein, " 2>/dev/null"))
      ## TODO: Take into account that CDO renames 'time' dimension to 'leadtime' or 'reftime' if these variables exist, and
      ##       probably the time axes get messed up. This is bad for file per dataset datasets.
      ## TODO: When loading file per dataset datasets and the leadtime dimension was removed by cdo (after bicubic interp.)
      ##       ncks -A -v leadtime /cfunas/exp/ENSEMBLES/decadal/MME/cerfacs/2002/s000/m001/monthly_mean/hfssd/hfssd_19601101.nc
    } else {
      filein <- filename
    }
    fnc <- nc_open(filein)
    var_dimnames <- names(fnc$dim)[fnc$var[[namevar]][['dimids']] + 1]
    nmemb <- nltime <- NULL
    if (work_piece[['dimnames']][['members']] %in% var_dimnames) {
      nmemb <- fnc$dim[[work_piece[['dimnames']][['members']]]]$len
    } else {
      nmemb <- 1
    }
    ## The steps below have to be done since CDO renames ''randomly'' the time dimension after remap.
    dim_matches <- match(c(work_piece[['dimnames']][['longitudes']],
                           work_piece[['dimnames']][['latitudes']],
                           work_piece[['dimnames']][['members']]), var_dimnames)
    dim_matches <- dim_matches[!is.na(dim_matches)]
    time_dimname <- var_dimnames[-dim_matches]
    if (length(time_dimname) > 0) {
      nltime <- fnc$dim[[time_dimname]]$len
    } else {
      nltime <- 1
    }
    
    if (work_piece[['is_2d_var']]) {
      ## We read the longitudes and latitudes from the file.
      lon <- ncvar_get(fnc, work_piece[['dimnames']][['longitudes']])
      full_grid_n_lons <- length(lon)
      lonmin <- work_piece[['lon_limits']][1]
      lonmax <- work_piece[['lon_limits']][2]
      lon[which(lon < 0)] <- lon[which(lon < 0)] + 360
      lon_indices <- 1:length(lon)
      lon_indices <- lon_indices[which((lon >= lonmin & lon <= lonmax) | 
                     (lonmax < lonmin & (lon >= lonmin | lon <= lonmax)))]
      lon <- lon[lon_indices]
      lat <- ncvar_get(fnc, work_piece[['dimnames']][['latitudes']])
      full_grid_n_lats <- length(lat)
      latmin <- work_piece[['lat_limits']][1]
      latmax <- work_piece[['lat_limits']][2]
      lat_indices <- 1:length(lat)
      if (lat[1] < lat[length(lat)]) {
        lat_indices <- length(lat):1
        lat <- lat[lat_indices]
      }
      lat_indices <- lat_indices[which(lat >= latmin & lat <= latmax)]
      lat <- lat[which(lat >= latmin & lat <= latmax)]
    } else {
      lon <- 0
      lat <- 0
    }

    if (explore_dims) {
      if (work_piece[['is_file_per_member']]) {
        ## TODO: When the exp_full_path contains asterisks and is file_per_member
        ##       members from different datasets may be accounted.
        ##       Also if one file member is missing the accounting will be wrong.
        ##       Should parse the file name and extract number of members.
        nmemb <- length(files)
      }
      dims <- list(nmemb = nmemb, nltime = nltime, lon = lon, lat = lat)
    } else {
      if (!is.null(work_piece[['mask']]) && !identical(c(full_grid_n_lons, full_grid_n_lats), dim(work_piece[['mask']]))) {
        stop(paste("Error: the mask of the dataset with index ", tail(work_piece[['indices']], 1), " in '", work_piece[['dataset_type']], "' is wrong. It must be on the common grid if the selected output type is 'lonlat', 'lon' or 'lat', or 'areave' and 'grid' has been specified. It must be on the grid of the corresponding dataset if the selected output type is 'areave' and no 'grid' has been specified. For more information check ?Load and see help on parameters 'grid', 'maskmod' and 'maskobs'.", sep = ""))
      }
      # The data is arranged in the array 'tmp' with the dimensions in a 
      # common order:
      #   1) Longitudes 
      #   2) Latitudes
      #   3) Members (even if is not a file per member experiment)
      #   4) Lead-times
      if (work_piece[['is_file_per_dataset']]) {
        time_indices <- 1:nltime
        mons <- strsplit(system(paste('cdo showmon ', filein, 
                         ' 2>/dev/null'), intern = TRUE), split = ' ')
        years <- strsplit(system(paste('cdo showyear ', filein, 
                          ' 2>/dev/null'), intern = TRUE), split = ' ')
        mons <- as.integer(mons[[1]][which(mons[[1]] != "")])
        years <- as.integer(years[[1]][which(years[[1]] != "")])
        time_indices <- ts(time_indices, start = c(years[1], mons[1]), 
                           end = c(years[length(years)], mons[length(mons)]),
                           frequency = 12)
        ltimes_list <- list()
        for (sdate in work_piece[['startdates']]) {
          selected_time_indices <- window(time_indices, start = c(as.integer(
                                   substr(sdate, 1, 4)), as.integer(substr(sdate, 5, 6))), 
                                   end = c(3000, 12), frequency = 12, extend = TRUE)
          ltimes_list <- c(ltimes_list, list(selected_time_indices[work_piece[['leadtimes']]]))
        }
      } else {
        ltimes <- work_piece[['leadtimes']]
        #if (work_piece[['dataset_type']] == 'exp') {
          ltimes_list <- list(ltimes[which(ltimes <= nltime)])
        #}
      }
      ## TODO: Put when reading matrices this kind of warnings
      #  if (nmember < nmemb) {
      #    cat("Warning:
      members <- 1:work_piece[['nmember']]
      members <- members[which(members <= nmemb)]

      for (ltimes in ltimes_list) {
        ## TODO: Pick only the needed data using ncdf4 extra parameters
        if (work_piece[['is_2d_var']]) {
          if (nmemb == 1 && nltime == 1) {
            tmp <- ncvar_get(fnc, namevar)[lon_indices, lat_indices]
            dim(tmp) <- c(length(lon_indices), length(lat_indices), 1, 1)
          } else if (nmemb == 1) {
            tmp <- ncvar_get(fnc, namevar)[lon_indices, lat_indices, ltimes]
            dim(tmp) <- c(length(lon_indices), length(lat_indices), 1, length(ltimes)) 
          } else if (nltime == 1) {
            tmp <- ncvar_get(fnc, namevar)[lon_indices, lat_indices, members]
            dim(tmp) <- c(length(lon_indices), length(lat_indices), length(members), 1)
          } else {
            tmp <- ncvar_get(fnc, namevar)[lon_indices, lat_indices, members, ltimes]
            dim(tmp) <- c(length(lon_indices), length(lat_indices), length(members), length(ltimes))
          }
          original_dims <- dim(tmp)
 
          mask <- work_piece[['mask']][lon_indices, lat_indices]
          tmp <- apply(tmp, c(3, 4), function(x) {
            # Disable of large values.
            if (!is.na(work_piece[['var_limits']][2])) {
              x[which(x > work_piece[['var_limits']][2])] <- NA
            }
            if (!is.na(work_piece[['var_limits']][1])) {
              x[which(x < work_piece[['var_limits']][1])] <- NA
            }
            if (!is.null(mask)) {
              x[which(mask < 0.5)] <- NA
            }

            if (output == 'areave' || output == 'lon') {
              weights <- InsertDim(cos(lat * pi / 180), 1, length(lon))
              weights[which(is.na(x))] <- NA
              if (output == 'areave') {
                weights <- weights / mean(weights, na.rm = TRUE)
                mean(x * weights, na.rm = TRUE) 
              } else {
                weights <- weights / InsertDim(Mean1Dim(weights, 2, narm = TRUE), 2, length(lat))
                Mean1Dim(x * weights, 2, narm = TRUE)
              }
            } else if (output == 'lat') {
              Mean1Dim(x, 1, narm = TRUE)
            } else if (output == 'lonlat') {
              signif(x, 5)
            }
          })
          if (output == 'areave') {
            dim(tmp) <- c(1, 1, original_dims[3:4])
          } else if (output == 'lon') {
            dim(tmp) <- c(original_dims[1], 1, original_dims[3:4])
          } else if (output == 'lat') {
            dim(tmp) <- c(1, original_dims[c(2, 3, 4)])
          } else if (output == 'lonlat') {
            dim(tmp) <- original_dims
          }
        } else {
          if (nmemb == 1 && nltime == 1) {
            tmp <- ncvar_get(fnc, namevar)
            dim(tmp) <- c(1, 1, 1, 1) 
          } else if (nmemb == 1) {
            tmp <- ncvar_get(fnc, namevar)[ltimes]
            dim(tmp) <- c(1, 1, 1, dim(tmp)) 
          } else if (nltime == 1) {
            tmp <- ncvar_get(fnc, namevar)[members]
            dim(tmp) <- c(1, 1, dim(tmp), 1) 
          } else {
            tmp <- ncvar_get(fnc, namevar)[members, ltimes]
            dim(tmp) <- c(1, 1, dim(tmp)) 
          }
        }
      
        var_data <- attach.big.matrix(work_piece[['out_pointer']])
        if (work_piece[['dims']][['nmember']] > 1 && nmemb > 1 && work_piece[['dims']][['nltime']] > 1 && nltime < work_piece[['dims']][['nltime']]) {
          work_piece[['indices']][2] <- work_piece[['indices']][2] - 1
          for (jmemb in 1:nmemb) {
            work_piece[['indices']][2] <- work_piece[['indices']][2] + 1
            out_position <- arrayIndex2VectorIndex(work_piece[['indices']], work_piece[['dims']])
            out_indices <- out_position:(out_position + length(tmp[, , jmemb, ]) - 1)
            var_data[out_indices] <- as.vector(tmp[, , jmemb, ])
          }
          work_piece[['indices']][2] <- work_piece[['indices']][2] - nmemb + 1
        } else {
          out_position <- arrayIndex2VectorIndex(work_piece[['indices']], work_piece[['dims']])
          out_indices <- out_position:(out_position + length(tmp) - 1)
          var_data[out_indices] <- as.vector(aperm(tmp, c(1, 2, 4, 3)))
        }

        work_piece[['indices']][3] <- work_piece[['indices']][3] + 1
      }
    }

    nc_close(fnc)    
    if (work_piece[['is_2d_var']] && !is.null(work_piece[['grid']])) {
      system(paste("rm -f", filein))
    }
  }

  if (explore_dims) {
    dims
  } else {
    file_found
  }
}

.LoadSampleData <- function(var, exp = NULL, obs = NULL, sdates, 
                            nmember = NULL, nmemberobs = NULL, 
                            nleadtime = NULL, leadtimemin = 1, 
                            leadtimemax = NULL, storefreq = 'monthly', 
                            sampleperiod = 1, lonmin = 0, lonmax = 360, 
                            latmin = -90, latmax = 90, output = 'areave', 
                            method = 'conservative', grid = NULL, 
                            maskmod = vector("list", 15), 
                            maskobs = vector("list", 15), 
                            configfile = NULL, suffixexp = NULL, 
                            suffixobs = NULL, varmin = NULL, varmax = NULL, 
                            silent = FALSE, nprocs = NULL) {
  ## This function loads and selects sample data stored in sampleMap and 
  ## sampleTimeSeries and is used in the examples instead of Load() so as
  ## to avoid nco and cdo system calls and computation time in the stage 
  ## of running examples in the CHECK process on CRAN.
  selected_start_dates <- match(sdates, c('19851101', '19901101', '19951101', 
                                          '20001101', '20051101'))
  start_dates_position <- 3
  lead_times_position <- 4

  if (output == 'lonlat') {
    sampleData <- s2dverification::sampleMap
    if (is.null(leadtimemax)) {
      leadtimemax <- dim(sampleData$mod)[lead_times_position]
    }
    selected_lead_times <- leadtimemin:leadtimemax

    dataOut <- sampleData
    dataOut$mod <- sampleData$mod[, , selected_start_dates, selected_lead_times, , ]
    dataOut$obs <- sampleData$obs[, , selected_start_dates, selected_lead_times, , ]
  }
  else if (output == 'areave') {
    sampleData <- s2dverification::sampleTimeSeries
    if (is.null(leadtimemax)) {
      leadtimemax <- dim(sampleData$mod)[lead_times_position]
    }
    selected_lead_times <- leadtimemin:leadtimemax

    dataOut <- sampleData
    dataOut$mod <- sampleData$mod[, , selected_start_dates, selected_lead_times]
    dataOut$obs <- sampleData$obs[, , selected_start_dates, selected_lead_times]
  }

  dims_out <- dim(sampleData$mod)
  dims_out[start_dates_position] <- length(selected_start_dates)
  dims_out[lead_times_position] <- length(selected_lead_times)
  dim(dataOut$mod) <- dims_out

  dims_out <- dim(sampleData$obs)
  dims_out[start_dates_position] <- length(selected_start_dates)
  dims_out[lead_times_position] <- length(selected_lead_times)
  dim(dataOut$obs) <- dims_out

  invisible(list(mod = dataOut$mod, obs = dataOut$obs, 
                 lat = dataOut$lat, lon = dataOut$lon))
}

.ConfigGetDatasetInfo <- function(matching_entries) {
  # This function obtains the information of a dataset and variable pair,
  # applying all the entries that match in the configuration file.
  experiment_defaults <- c('$DEFAULT_EXP_MAIN_PATH$', '$DEFAULT_EXP_FILE_PATH$', '$DEFAULT_GRID$', '$DEFAULT_NC_VAR_NAME$', '$DEFAULT_SUFFIX$', '$DEFAULT_VAR_MIN$', '$DEFAULT_VAR_MAX$')
  observation_defaults <- c('$DEFAULT_OBS_MAIN_PATH$', '$DEFAULT_OBS_FILE_PATH$', '$DEFAULT_NC_VAR_NAME$', '$DEFAULT_SUFFIX$', '$DEFAULT_VAR_MIN$', '$DEFAULT_VAR_MAX$')
  info <- NULL

  for (entry in matching_entries) {
    if (is.null(info)) {
      info <- entry[-1:-2]
      if (length(info) == 7) {
        info[which(info == '*')] <- experiment_defaults[which(info == '*')]
      } else {
        info[which(info == '*')] <- observation_defaults[which(info == '*')]
      }
    } else {
      info[which(entry[-1:-2] != '*')] <- entry[-1:-2][which(entry[-1:-2] != '*')]
    }
  }
  
  info
}
