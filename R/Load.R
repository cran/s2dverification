Load <- function(var, exp = NULL, obs = NULL, sdates, nmember = NULL, 
                 nmemberobs = NULL, nleadtime = NULL, leadtimemin = 1, 
                 leadtimemax = NULL, storefreq = 'monthly', sampleperiod = 1, 
                 lonmin = 0, lonmax = 360, latmin = -90, latmax = 90, 
                 output = 'areave', method = 'conservative', grid = NULL, 
                 maskmod = vector("list", 15), maskobs = vector("list", 15), 
                 configfile = NULL, suffixexp = NULL, suffixobs = NULL, 
                 varmin = NULL, varmax = NULL, silent = FALSE, nprocs = NULL,
                 dimnames = NULL) {
  #library(parallel)
  #library(bigmemory)

  # Check parameters
  if (!(storefreq %in% c('monthly', 'daily'))) {
    stop("Error: parameter 'storefreq' is wrong, can take value 'daily' or 'monthly'.")
  }
  if (!(method %in% c('bilinear', 'bicubic', 'conservative', 'distance-weighted'))) {
    stop("Error: parameter 'method' is wrong, can take value 'bilinear', 'bicubic', 'conservative' or 'distance-weighted'.")
  }
  if (is.null(configfile)) {
    configfile <- system.file("config", "IC3.conf", package = "s2dverification")
  }
  if (!is.null(suffixexp)) {
    if (length(suffixexp) == 1) {
      suffixexp <- rep(suffixexp[1], length(exp))
    }
    if (length(suffixexp) != length(exp)) {
      stop("Error: parameter 'suffixexp' must contain 0, 1 or length(exp) values.")
    }
  }
  if (!is.null(suffixobs)) {
    if (length(suffixobs) == 1) {
      suffixobs <- rep(suffixobs[1], length(obs))
    }
    if (length(suffixobs) != length(obs)) {
      stop("Error: parameter 'suffixobs' must contain 0, 1 or length(obs) values.")
    }
  }
  if (!is.null(nmember) && !is.null(exp)) {
    if (length(nmember) == 1) {
      cat(paste("Warning: 'nmember' should specify the number of members of each experimental dataset. Forcing to", nmember, "for all experiments.\n"))
      nmember <- rep(nmember, length(exp))
    }
    if (length(nmember) != length(exp)) {
      stop("Error: 'nmember' must contain as many values as 'exp'.")
    } else if (any(is.na(nmember))) {
      nmember[which(is.na(nmember))] <- max(nmember, na.rm = TRUE)
    }
  }
  if (!is.null(nmemberobs) && !is.null(obs)) {
    if (length(nmemberobs) == 1) {
      cat(paste("Warning: 'nmemberobs' should specify the number of members of each observational dataset. Forcing to", nmemberobs, "for all observations.\n"))
      nmemberobs <- rep(nmemberobs, length(obs))
    }
    if (length(nmemberobs) != length(obs)) {
      stop("Error: 'nmemberobs' must contain as many values as 'obs'.")
    } else if (any(is.na(nmemberobs))) {
      nmemberobs[which(is.na(nmemberobs))] <- max(nmemberobs, na.rm = TRUE)
    }
  }
  remap <- switch(method, 'bilinear' = 'remapbil', 'bicubic' = 'remapbic', 
                  'conservative' = 'remapcon', 'distance-weighted' = 'remapdis')
  if (!(output %in% c('lonlat', 'lon', 'lat', 'areave'))) {
    stop(paste("Error: 'output' can only take values 'lonlat', 'lon', 'lat' or 'areave'."))
  }
  ## Make sure that longitude bounds are positive
  if (lonmin < 0) {
    lonmin <- lonmin + 360
  }
  if (lonmax < 0) {
    lonmax <- lonmax + 360
  }
  ## Force the observational masks to be the same as the experimental when
  ## possible.
  if (length(maskmod) < length(exp)) {
    stop("Error: 'maskmod' must contain a mask for each experiment in 'exp'.")
  }
  if ((output != 'areave' || !is.null(grid)) && length(exp) > 0) {
    if (!identical(maskobs, vector("list", 15))) {
      cat("Warning: 'maskobs' will be ignored. 'maskmod[[1]]' will be applied to observations instead\n.")
    }
    maskobs <- lapply(maskobs, function(x) x <- maskmod[[1]])
  } else if (length(maskobs) < length(obs)) {
    stop("Error: 'maskobs' must contain a mask for each observation in 'obs'.")
  }

  # Extend the current environment with the functions and variables related to 
  # the configuration file mechanism which contain the information about the 
  # location of the datasets.
  # 'replace_values' is a named list that associates a variable name to an 
  # associated value. Initially it is filled with variables and values parsed
  # from the configuration file, but we can add or modify some values during
  # the execution to choose for example which start date we want to load.
  # When '.ConfigReplaceVariablesInString' is called, all the variable accesses 
  # ($VARIABLE_NAME$) that appear in the string given as parameter are 
  # replaced by the associated value in 'replace_values'.
  #
  # ConfigFileOpen will tell if the variable we are loading is 
  # 2-dimensional or not and if the datasets we are loading are stored in
  # file per member, file per ensemble or file per dataset format.
  data_info <- ConfigFileOpen(configfile, silent)
  
  # Check that the var, exp and obs parameters are right and keep the entries that match for each dataset.
  # Afterwards, the matching entries are applied sequentially (as specified in ?ConfigFileOpen)
  # and the replace_values are applied to the result.
  # Finally a path pattern for each dataset is provided.
  is_2d_var <- ConfigMatchingVars(data_info, var, show = FALSE)
  matches <- ConfigApplyMatchingEntries(data_info, var, exp, obs, show_entries = FALSE, show_result = FALSE)
  exp_info <- matches$exp_info
  is_file_per_member_exp <- matches$is_file_per_member_exp
  obs_info <- matches$obs_info
  is_file_per_member_obs <- matches$is_file_per_member_obs
  is_file_per_dataset_obs <- matches$is_file_per_dataset_obs
  
  replace_values <- data_info$definitions
  ## We add in this table of variables the predefined values.
  replace_values[["VAR_NAME"]] <- var
  replace_values[["STORE_FREQ"]] <- storefreq

  if (!silent) {
    cat("* All pairs (var, exp) and (var, obs) have matching entries.\n")
  }

  # Some extra checks now that we know the variable type
  if (!is_2d_var && (output != 'areave')) {
    cat(paste("Warning: '", output, "' output format not allowed when loading global mean variables. Forcing to 'areave'\n.", 
                sep = ''))
    output <- 'areave'
  }

  # Work out the nc file dimension names
  dimnames <- list(longitudes = ifelse(is.null(dimnames[["longitudes"]]), 
                                       replace_values[["DEFAULT_DIM_NAME_LONGITUDES"]],
                                       dimnames[['longitudes']]),
                   latitudes = ifelse(is.null(dimnames[["latitudes"]]), 
                                       replace_values[["DEFAULT_DIM_NAME_LATITUDES"]],
                                       dimnames[['latitudes']]),
                   members = ifelse(is.null(dimnames[["members"]]), 
                                       replace_values[["DEFAULT_DIM_NAME_MEMBERS"]],
                                       dimnames[['members']]))

  # Initialize some variables that will take various values along the
  # execution
  latitudes <- longitudes <- NULL
  leadtimes <- NULL
  var_exp <- var_obs <- NULL

  # Start defining the dimensions of the output matrices
  nmod <- length(exp)
  nobs <- length(obs)
  nsdates <- length(sdates)

  # We will iterate over all the experiments, start dates and members and will open
  # the file pointed by the data in the configuration file.
  # If a file is found, we will open it and read its metadata to work out the 
  # remaining dimensions: members, leadtimes, longitudes and latitudes.
  #
  # At each iteration we will build a 'work piece' that will contain information
  # on the data we want to load from a file. For each file we will have one
  # work piece. These work pieces will be packages of information to be sent to
  # the various parallel processes. Each process will need this information to
  # access and manipulate the data according to the output type and other 
  # parameters. 
  if (!silent) {
    cat("* Fetching first experimental files to work out 'var_exp' size...\n")
  }

  dataset_type <- 'exp'
  dim_exp <- NULL
  filename <- file_found <- tmp <- nltime <- NULL
  dims2define <- TRUE
  exp_work_pieces <- list()
  jmod <- 1
  while (jmod <= nmod) {
    replace_values[["EXP_NAME"]] <- exp[jmod]
    dataset_info <- exp_info[[jmod]]
    replace_values[["EXP_MAIN_PATH"]] <- dataset_info[1]
    replace_values[["EXP_FILE_PATH"]] <- dataset_info[2]
    replace_values[["GRID"]] <- dataset_info[3]
    replace_values[["NC_VAR_NAME"]] <- dataset_info[4]
    namevar <- .ConfigReplaceVariablesInString(dataset_info[4], replace_values)
    if (is.null(suffixexp) || is.na(suffixexp[jmod])) {
      replace_values[["SUFFIX"]] <- dataset_info[5]
    } else {
      replace_values[["SUFFIX"]] <- suffixexp[jmod]
    }
    if (is.null(varmin)) {
      mod_var_min <- as.numeric(.ConfigReplaceVariablesInString(dataset_info[6], replace_values))
    } else {
      mod_var_min <- varmin
    }
    if (is.null(varmax)) {
      mod_var_max <- as.numeric(.ConfigReplaceVariablesInString(dataset_info[7], replace_values))
    } else {
      mod_var_max <- varmax
    }
    jsdate <- 1
    while (jsdate <= nsdates) {
      replace_values[["START_DATE"]] <- sdates[jsdate]
      replace_values[["YEAR"]] <- substr(jsdate, 1, 4)
      replace_values[["MONTH"]] <- substr(jsdate, 5, 6)
      replace_values[["DAY"]] <- substr(jsdate, 7, 8)
      # If the dimensions of the output matrices are still to define, we try to read
      # the metadata of the data file that corresponds to the current iteration
      if (dims2define) {
        if (is_file_per_member_exp[jmod]) {
          replace_values[["MEMBER_NUMBER"]] <- '*'
        }
        # Now we check that the specified grid is supported by CDO
        grid_to_try <- grid
        if (output != 'areave') {
          if (is.null(grid)) {
            grid_to_try <- .ConfigReplaceVariablesInString(dataset_info[3], replace_values)
          }
          supported_grids <- list('r[0-9]{1,}x[0-9]{1,}', 't[0-9]{1,}grid')
          error_message <- paste('The specified grid in', ifelse(is.null(grid), paste0('the configuration file for the experiment ', exp[jmod]), "the parameter 'grid'"), 'is incorrect. Must be one of r<NX>x<NY> or t<RES>grid.')
          grid_matches <- unlist(lapply(lapply(supported_grids, regexpr, grid_to_try), .IsFullMatch, grid_to_try))
          if (sum(grid_matches) < 1) {
            stop(error_message)
          }
          ## TODO: If is right, we should create and pass an object with grid info in work_piece,
          ##       in place of 'grid'.
          ##       This way it will be possible to check if interpolation is unnecessary.
        }
        # We must build a work piece that will be sent to the .LoadDataFile function
        # in 'explore_dims' mode. We will obtain, if success, the dimensions of the
        # data in the file.
        work_piece <- list(dataset_type = dataset_type,
                           filename = .ConfigReplaceVariablesInString("$EXP_FULL_PATH$", replace_values),
                           namevar = namevar, is_2d_var = is_2d_var, grid = grid_to_try,
                           remap = remap,
                           is_file_per_member = is_file_per_member_exp[jmod],
                           lon_limits = c(lonmin, lonmax),
                           lat_limits = c(latmin, latmax), dimnames = dimnames)
        found_dims <- .LoadDataFile(work_piece, explore_dims = TRUE, silent = silent)
        if (!is.null(found_dims)) {
          grid <- grid_to_try
          if (is.null(nmember)) {
            nmember <- rep(found_dims[['nmemb']], nmod)
          }
          if (is.null(nleadtime)) {
            nleadtime <- found_dims[['nltime']]
          }
          if (is.null(leadtimemax)) {
            leadtimemax <- nleadtime
          } else if (leadtimemax > nleadtime) {
            stop("Error: 'leadtimemax' argument is greater than the number of loadedleadtimes. 
                  Put first the experiment with the greatest number of leadtimes or adjust the parameters 'nleadtime' and 'leadtimemax' properly.")
          }
          leadtimes <- seq(leadtimemin, leadtimemax, sampleperiod)
          latitudes <- found_dims[['lat']]
          longitudes <- found_dims[['lon']]
          
          if (output == 'lon' || output == 'lonlat') {
            dim_exp[['lon']] <- length(longitudes)
          }
          if (output == 'lat' || output == 'lonlat') {
            dim_exp[['lat']] <- length(latitudes)
          }
          dim_exp[['nltime']] <- length(leadtimes)
          dim_exp[['nmember']] <- max(nmember)
          dim_exp[['nsdates']] <- nsdates
          dim_exp[['ndat']] <- nmod
          dims2define <- FALSE
        }
      }
      # We keep on iterating through members to build all the work pieces.
      if (is_file_per_member_exp[jmod]) {
        jmember <- 1
        while (jmember <= nmember[jmod]) {
          replace_values[["MEMBER_NUMBER"]] <- sprintf(paste("%.", (nmember[jmod] %/% 10) + 1, "i", sep = ''), jmember - 1)
          work_piece <- list(filename = .ConfigReplaceVariablesInString("$EXP_FULL_PATH$", replace_values),
                             namevar = namevar, indices = c(1, jmember, jsdate, jmod), 
                             nmember = nmember[jmod], leadtimes = leadtimes, mask = maskmod[[jmod]],
                             is_file_per_dataset = FALSE,
                             var_limits = c(mod_var_min, mod_var_max))
          exp_work_pieces <- c(exp_work_pieces, list(work_piece))
          jmember <- jmember + 1
        }
      } else {
        work_piece <- list(filename = .ConfigReplaceVariablesInString("$EXP_FULL_PATH$", replace_values),
                           namevar = namevar, indices = c(1, 1, jsdate, jmod), 
                           nmember = nmember[jmod], leadtimes = leadtimes, mask = maskmod[[jmod]],
                           is_file_per_dataset = FALSE,
                           var_limits = c(mod_var_min, mod_var_max))
        exp_work_pieces <- c(exp_work_pieces, list(work_piece))
      }
      jsdate <- jsdate + 1
    }
    jmod <- jmod + 1
  }
  if (dims2define && length(exp) > 0) {
    cat("Warning: no data found in file system for any experimental dataset. Forcing 'exp' to NULL.\n")
    exp <- NULL
    dims2define <- FALSE
  }
  # If there are no experiments to load we need to choose a number of time steps
  # to load from observational datasets. We load from the first start date to 
  # the current date.
  if (is.null(exp)) {
    diff <- Sys.time() - as.POSIXct(paste(substr(sdates[1], 1, 4), '-',
            substr(sdates[1], 5, 6), '-', substr(sdates[1], 7, 8), sep=''))
    if (storefreq == 'monthly') { 
      leadtimemax <- as.integer(diff/30) 
    } else {
      leadtimemax <- as.integer(diff)
    }
    if (is.null(nleadtime)) {
      nleadtime <- leadtimemax
    }
    leadtimes <- seq(leadtimemin, leadtimemax, sampleperiod)
  }

  if (!silent) {
    message <- paste("* Success: ", nmod, ", ", 
                     ifelse(is.null(nmember), 0, max(nmember)), ", ", 
                     nsdates, ", ", length(leadtimes), sep = '')
    if (output == 'lon') {
      message <- paste(message, ", ", length(longitudes), sep = '')
    } else if (output == 'lat') {
      message <- paste(message, ", ", length(longitudes), sep = '')
    } else if (output == 'lonlat') {
      message <- paste(message, ", ", length(latitudes), ", ", length(longitudes), sep = '')
    }
    cat(paste0(message ,'\n'))
    cat("* Fetching first observational files to work out 'var_obs' size...\n")
  }
  
  # Now we start iterating over observations. We try to find the output matrix
  # dimensions and we build anyway the work pieces corresponding to the observational
  # data that time-corresponds the experimental data or the time-steps until the
  # current date if no experimental datasets were specified.
  dataset_type <- 'obs'
  dim_obs <- NULL
  dims2define <- TRUE
  lat_indices <- lon_indices <- NULL
  obs_work_pieces <- list()
  jobs <- 1
  while (jobs <= nobs) {
    replace_values[["OBS_NAME"]] <- obs[jobs]
    dataset_info <- obs_info[[jobs]]
    replace_values[["OBS_MAIN_PATH"]] <- dataset_info[1]
    replace_values[["OBS_FILE_PATH"]] <- dataset_info[2]
    replace_values[["NC_VAR_NAME"]] <- dataset_info[3]
    namevar <- .ConfigReplaceVariablesInString(dataset_info[3], replace_values)
    if (is.null(suffixobs) || is.na(suffixobs[jobs])) {
      replace_values[["SUFFIX"]] <- dataset_info[4]
    } else {
      replace_values[["SUFFIX"]] <- suffixobs[jobs]
    }
    if (is.null(varmin)) {
      obs_var_min <- as.numeric(.ConfigReplaceVariablesInString(dataset_info[5], replace_values))
    } else {
      obs_var_min <- varmin
    }
    if (is.null(varmax)) {
      obs_var_max <- as.numeric(.ConfigReplaceVariablesInString(dataset_info[6], replace_values))
    } else {
      obs_var_max <- varmax
    }
    # This file format (file per whole dataset) is only supported in observations.
    # However a file per whole dataset experiment could be seen as a file per
    # member/ensemble experiment with a single start date, so still loadable.
    # Nonetheless file per whole dataset observational files do not need to contain
    # a year and month in the filename, the time correspondance relies on the 
    # month and years associated to each timestep inside the NetCDF file.
    # So file per whole dataset experiments need to have a start date in the filename.
    if (is_file_per_dataset_obs[jobs]) {
      ## TODO: Open file-per-dataset-files only once.
      if (dims2define) {
        work_piece <- list(dataset_type = dataset_type,
                           filename = .ConfigReplaceVariablesInString("$OBS_FULL_PATH$", replace_values),
                           namevar = namevar, is_2d_var = is_2d_var, grid = grid,
                           remap = remap, is_file_per_member = is_file_per_member_obs[jobs],
                           lon_limits = c(lonmin, lonmax),
                           lat_limits = c(latmin, latmax), dimnames = dimnames)
        found_dims <- .LoadDataFile(work_piece, explore_dims = TRUE, silent = silent)
        if (!is.null(found_dims)) {
          if (is.null(nmemberobs)) {
            nmemberobs <- rep(found_dims[['nmemb']], nobs)
          }
          if (is.null(exp)) {
            longitudes <- found_dims[['lon']]
            latitudes <- found_dims[['lat']]
          }
          
          if (output == 'lon' || output == 'lonlat') {
            dim_obs[['lon']] <- length(longitudes)
          }
          if (output == 'lat' || output == 'lonlat') {
            dim_obs[['lat']] <- length(latitudes)
          }
          dim_obs[['nltime']] <- length(leadtimes)
          dim_obs[['nmember']] <- max(nmemberobs)
          dim_obs[['nsdates']] <- nsdates
          dim_obs[['ndat']] <- nobs
          dims2define <- FALSE
        }
      }
    
      work_piece <- list(filename = .ConfigReplaceVariablesInString("$OBS_FULL_PATH$", replace_values),
                         namevar = namevar, indices = c(1, 1, 1, jobs), 
                         nmember = nmemberobs[jobs], 
                         mask = maskobs[[jobs]], leadtimes = leadtimes, 
                         is_file_per_dataset = is_file_per_dataset_obs[jobs], 
                         startdates = sdates, var_limits = c(obs_var_min, obs_var_max))
      obs_work_pieces <- c(obs_work_pieces, list(work_piece))
    } else {
      jsdate <- 1
      while (jsdate <= nsdates) {
        replace_values[["START_DATE"]] <- sdates[jsdate]        
        sdate <- sdates[jsdate]

        if (storefreq == 'daily') {
          day <- substr(sdate, 7, 8)
          if (day == '') {
            day <- '01'
          }
          day <- as.integer(day)
          startdate <- as.POSIXct(paste(substr(sdate, 1, 4), '-',
                       substr(sdate, 5, 6), '-', day, ' 12:00:00', sep = '')) + 
                       (leadtimemin - 1) * 86400
          year <- as.integer(substr(startdate, 1, 4))
          month <- as.integer(substr(startdate, 6, 7))
        } else {
          month <- (as.integer(substr(sdate, 5, 6)) + leadtimemin - 2) %% 12 + 1
          year <- as.integer(substr(sdate, 1, 4)) + (as.integer(substr(sdate, 
                  5, 6)) + leadtimemin - 2) %/% 12
        }
        jleadtime <- 1
        while (jleadtime <= length(leadtimes)) {
          replace_values[["YEAR"]] <- paste(year, '', sep = '')
          replace_values[["MONTH"]] <- sprintf("%2.2i", month)
        
          if (storefreq == 'daily') {
            replace_values[["DAY"]] <- sprintf("%2.2i", day)
            days_in_month <- ifelse(LeapYear(year), 29, 28)
            days_in_month <- switch(paste(month, '', sep = ''), '1' = 31, 
                                    '3' = 31, '4' = 30, '5' = 31, '6' = 30, 
                                    '7' = 31, '8' = 31, '9' = 30, '10' = 31, 
                                    '11' = 30, '12' = 31, days_in_month)
            ## This condition must be fulfilled to put all the month time steps
            ## in the dimension of length nleadtimes. Otherwise it must be cut:
            #(length(leadtimes) - 1)*sampleperiod + 1 - (jleadtime - 1)*sampleperiod >= days_in_month - day + 1
            obs_file_indices <- seq(day, min(days_in_month, (length(leadtimes) - jleadtime) * sampleperiod + day), sampleperiod)
          } else {
            obs_file_indices <- 1
          }
          if (dims2define) {
            if (is_file_per_member_obs[jobs]) {
              replace_values[["MEMBER_NUMBER"]] <- '*'
            }
            work_piece <- list(dataset_type = dataset_type,
                               filename = .ConfigReplaceVariablesInString("$OBS_FULL_PATH$", replace_values),
                               namevar = namevar, is_2d_var = is_2d_var, grid = grid,
                               remap = remap,
                               is_file_per_member = is_file_per_member_obs[jobs],
                               lon_limits = c(lonmin, lonmax),
                               lat_limits = c(latmin, latmax), dimnames = dimnames)
            found_dims <- .LoadDataFile(work_piece, explore_dims = TRUE, silent = silent)
            if (!is.null(found_dims)) {
              if (is.null(nmemberobs)) {
               nmemberobs <- rep(found_dims[['nmemb']], nobs)
              }
              if (is.null(exp)) {
                longitudes <- found_dims[['lon']]
                latitudes <- found_dims[['lat']]
              }
              
              if (output == 'lon' || output == 'lonlat') {
                dim_obs[['lon']] <- length(longitudes)
              }
              if (output == 'lat' || output == 'lonlat') {
                dim_obs[['lat']] <- length(latitudes)
              }
              dim_obs[['nltime']] <- length(leadtimes)
              dim_obs[['nmember']] <- max(nmemberobs)
              dim_obs[['nsdates']] <- nsdates
              dim_obs[['ndat']] <- nobs
              dims2define <- FALSE
            }
          }
          if (is_file_per_member_obs[jobs]) {
            jmember <- 1
            while (jmember <= nmemberobs[jobs]) {
              replace_values[["MEMBER_NUMBER"]] <- sprintf(paste("%.", (nmemberobs[jobs] %/% 10) + 1, "i", sep = ''), jmember - 1)
              work_piece <- list(filename = .ConfigReplaceVariablesInString("$OBS_FULL_PATH$", replace_values),
                                 namevar = namevar, indices = c(jleadtime, jmember, jsdate, jobs), 
                                 nmember = nmemberobs[jobs], leadtimes = obs_file_indices, 
                                 mask = maskobs[[jobs]],
                                 is_file_per_dataset = is_file_per_dataset_obs[jobs], 
                                 var_limits = c(obs_var_min, obs_var_max))
              obs_work_pieces <- c(obs_work_pieces, list(work_piece))
              jmember <- jmember + 1
            }
          } else {
            work_piece <- list(filename = .ConfigReplaceVariablesInString("$OBS_FULL_PATH$", replace_values),
                               namevar = namevar, indices = c(jleadtime, 1, jsdate, jobs), 
                               nmember = nmemberobs[jobs], leadtimes = obs_file_indices, 
                               mask = maskobs[[jobs]],
                               is_file_per_dataset = is_file_per_dataset_obs[jobs], 
                               var_limits = c(obs_var_min, obs_var_max))
            obs_work_pieces <- c(obs_work_pieces, list(work_piece))
          }
          
          if (storefreq == 'daily') {
            startdate <- startdate + 86400 * sampleperiod * length(obs_file_indices)
            year <- as.integer(substr(startdate, 1, 4))
            month <- as.integer(substr(startdate, 6, 7))
            day <- as.integer(substr(startdate, 9, 10))
          } else {
            month <- month + sampleperiod
            year <- year + (month - 1) %/% 12
            month <- (month - 1) %% 12 + 1
          }
          jleadtime <- jleadtime + length(obs_file_indices)
        }
        
        jsdate <- jsdate + 1
      }
    }
    jobs <- jobs + 1
  }
  if (dims2define && length(obs) > 0) {
    cat("Warning: no data found in file system for any observational dataset. Forcing 'obs' to NULL.\n")
    obs <- NULL
    dims2define <- FALSE
  }

  if (!silent) {
    message <- paste("* Success: ", nobs, ", ", 
                     ifelse(is.null(nmemberobs), 0, max(nmemberobs)), ", ", 
                     nsdates, ", ", length(leadtimes), sep = '')
    if (output == 'lon') {
      message <- paste(message, ", ", length(longitudes), sep = '')
    } else if (output == 'lat') {
      message <- paste(message, ", ", length(latitudes), sep = '')
    } else if (output == 'lonlat') {
      message <- paste(message, ", ", length(latitudes), ", ", 
                                      length(longitudes), sep = '')
    }
    cat(paste0(message, '\n'))
  }

  # We buil two matrices in shared memory for the parallel processes to
  # store their results
  # These matrices will contain data arranged with the following
  # dimension order, to maintain data spacial locality during the 
  # parallel fetch:
  #   longitudes, latitudes, leadtimes, members, startdates, nmod/nobs
  # So [1, 1, 1, 1, 1, 1] will be next to [2, 1, 1, 1, 1, 1] in memory
  pointer_var_exp <- NULL
  if (!is.null(exp)) {
    var_exp <- big.matrix(nrow = prod(dim_exp), ncol = 1)
    pointer_var_exp <- describe(var_exp)
  }
  pointer_var_obs <- NULL
  if (!is.null(obs)) {
    var_obs <- big.matrix(nrow = prod(dim_obs), ncol = 1)
    pointer_var_obs <- describe(var_obs)
  }

  if (is.null(nprocs)) {
    nprocs <- detectCores()
  }

  ## TODO: Try 2 cores only
  # Add some important extra fields in the work pieces before sending
  exp_work_pieces <- lapply(exp_work_pieces, function (x) c(x, list(dataset_type = 'exp', dims = dim_exp, out_pointer = pointer_var_exp)))
  obs_work_pieces <- lapply(obs_work_pieces, function (x) c(x, list(dataset_type = 'obs', dims = dim_obs, out_pointer = pointer_var_obs)))
  work_pieces <- c(exp_work_pieces, obs_work_pieces)
  if (!silent) {
    cat(paste("* Reading and processing data files (on ", nprocs, " parallel processes):\n", sep = ''))
  }
  work_pieces <- lapply(work_pieces, 
                 function (x) {
                   if (!silent) {
                     cat(paste("*   ", x[['filename']], '\n', sep = ''))
                   }
                   c(x, list(is_2d_var = is_2d_var, grid = grid, remap = remap,  
                             lon_limits = c(lonmin, lonmax), 
                             lat_limits = c(latmin, latmax), 
                             dimnames = dimnames, output = output))
                 })
  # Build the cluster of processes that will do the work and dispatch work pieces.
  # The function .LoadDataFile is applied to each work package. This function will
  # open the data file, regrid if needed, trim (select time steps, longitudes, 
  # latitudes, members), apply the mask, compute and apply the weights if needed,
  # disable extreme values and store in the shared memory matrices. 
  if (nprocs == 1) {
    found_files <- lapply(work_pieces, .LoadDataFile, silent = silent)
  } else {
    cluster <- makeCluster(nprocs, outfile = "")
    found_files <- clusterApply(cluster, work_pieces, .LoadDataFile, silent = silent)
    stopCluster(cluster)
  }
  if (!silent && any(!unlist(found_files))) {
    cat("* Warning: The following files were not found in the file system. Filling with NA values instead:\n")
    lapply(work_pieces[!unlist(found_files)], function (x) cat(paste("*   ", x[['filename']], '\n', sep = '')))
  }

  # Before ending, the data is arranged in the common format, with the following
  # dimension order:
  #  nmod/nobs, members, startdates, leadtimes, latitudes, longitudes
  invisible(list(mod = if (is.null(exp)) {
                         NULL
                       } else { 
                         dim_reorder <- length(dim_exp):1
                         dim_reorder[2:3] <- dim_reorder[3:2]
                         aperm(array(bigmemory::as.matrix(var_exp), dim = dim_exp), dim_reorder)
                       }, 
                 obs = if (is.null(obs)) {
                         NULL
                       } else {
                         dim_reorder <- length(dim_obs):1
                         dim_reorder[2:3] <- dim_reorder[3:2]
                         aperm(array(bigmemory::as.matrix(var_obs), dim = dim_obs), dim_reorder)
                       }, 
                 lat = latitudes, lon = longitudes))
}
