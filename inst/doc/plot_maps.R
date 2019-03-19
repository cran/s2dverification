#!/usr/bin/env Rscript

library(s2dverification)
library(ncdf4)
library(abind)
args <- commandArgs(TRUE)

comptrend <- TRUE    # Trend as a function of the start date for each leadtime
compcor <- TRUE      # Correlation Coefficient
comprms <- TRUE      # Root Mean Square Error
comprmsss <- TRUE    # Root Mean Square Skill Score
compratrms <- TRUE   # Ratio RMSE expid1 / expid2

var <- args[1]    # tos/tas/pr
#season <- args[2] # Year/DJF/MAM/JJA/SON
season <- 'all'
#ly <- args[3]     # ly1/ly2-5/ly6-9 for Forecast year 1 / years 2 to 5 / years 
                  # 6 to 9
ly <- 'all'
nltimemax <- 124  # number of leadtimes max in the experiments (in months)
nltimechunk <- 4  # number of leadtimes per chunk (in months)
lstexpid <- c('i00k','b02p') # list of ids
grid <- '320x160' # atmospheric grid for tos/pr (ocean/land only)
mon0 <- 11        # initial month
year0 <- 1960     # first start date
yearf <- 2005     # last start date
intsdate <- 5     # interval between start dates
runmeanlen <- 12  # length of the window for running mean (in months)
members <- list('19900101' = c('r1i1p1')) 
#PUT IN ORDER. NONE CAN BE MISSING
chunks <- list('19900101' = c('199001-199001', '199002-199002', '199003-199003'))



obs <- switch(var, 'tas' = 'ghcnersstgiss', 'tos' = 'ersst_v3b', 'pr' = 'gpcc1x1_v6')

syears <- seq(year0, yearf, intsdate)
imon2 <- paste("0", as.character(mon0), sep = "")
sdates <- paste(as.character(syears), substr(imon2, nchar(imon2) - 1, 
                nchar(imon2)), '01', sep = "")

savename <- paste(var, '_', season, '_', ly, sep = '')
for (expid in lstexpid ) {
  savename <- paste(savename, '_', expid, sep = '')
}

if (file.exists(paste(savename, '.sav', sep = ''))) {
  cat(paste0("Loading existing data from backup file ", savename, "...\n"))
  load(paste(savename, '.sav', sep = ''))
  if (is.null(toto1$mod) || is.null(toto1$obs)) {
    stop("Missing experimental or observational data in backup file. Remove the backup to try retrieving again.")
  }
} else {
  latbnd <- switch(var, 'tos' = c(-60, 65), 'pr' = c(-60, 65), c(-90, 90))
  varname <- switch(var, 'siarea' = paste(var, tolower(pole), sep = ''), 
                         'siextent' = paste(var, tolower(pole), sep = ''), 
                         'sivol' = paste(var, tolower(pole), sep = ''), 
                          var)
  #if (is.na(match('b02p', lstexpid)) == TRUE) {
    #lstload <- lstexpid
  #} else {
    #lstload <- lstexpid[-match('b02p', lstexpid)]
  #}
  exp <- list(path = paste0('/es*/exp/ecearth/', lstexpid[1], '/monthly_mean/$VAR_NAME$*/$VAR_NAME$_*_S$START_DATE$_$MEMBER$_$CHUNK$.nc'))
  b02p <- list(name = 'b02p')
  lstobs <- obs


  toto1 <- list()
  leading_empty_sdates <- 0
  greatest_n_of_members <- max(sapply(members, length))
  greatest_n_of_ltimes <- max(sapply(chunks, length)) * nltimechunk
  if (any(sapply(chunks, length) == 1)) {
	greatest_n_of_ltimes <- max(nltimeout, greatest_n_of_ltimes)
  }
  for (sdate in sdates) {
    sdate_number <- which(sdates == sdate)
	members_exist <- FALSE
	chunks_exist <- FALSE
	if (sdate %in% names(members)) {
      if (length(members[[sdate]]) > 0) {
		members_exist <- TRUE
	  }
	}
	if (sdate %in% names(chunks)) {
	  if (length(chunks[[sdate]]) > 0) {
		chunks_exist <- TRUE
	  }
	}
	load_exp_data <- FALSE
    toto1_exp_sdate <- NULL
	if (members_exist && chunks_exist) {
      load_exp_data <- TRUE
	} else {
      if (is.null(toto1$mod)) {
		leading_empty_sdates <- leading_empty_sdates + 1
	  } else {
		sdate_dims <- dim(toto1$mod)
		sdate_dims[3] <- 1
	    toto1_exp_sdate <- array(dim = sdate_dims)
      }
	}
	load_obs_data <- FALSE
	if (!is.null(lstobs)) {
	  load_obs_data <- TRUE
	}
    ltimes_to_take <- ifelse(length(chunks[[sdate]]) == 1, nltimeout, nltimechunk)
	toto1_obs_sdate <- NULL
    for (member in members[[sdate]]) {
      member_number <- which(members[[sdate]] == member)
      if (load_obs_data && member_number > 1) {
        load_obs_data <- FALSE
      }
      toto1_exp_member <- NULL
      for (chunk in chunks[[sdate]]) {
        chunk_number <- which(chunks[[sdate]] == chunk)
		if (load_exp_data) {
          lstload <- list(exp)
          lstload[[1]][['path']] <- gsub('$MEMBER$', member, lstload[[1]][['path']], fixed = TRUE)
          lstload[[1]][['path']] <- gsub('$CHUNK$', chunk, lstload[[1]][['path']], fixed = TRUE)
          toto1_exp_chunk <- Load(varname, lstload, NULL, sdate, leadtimemax = ltimes_to_take, dimnames = list(member = 'region'), output = 'lonlat', grid = grid)
          if (is.null(toto1_exp_chunk$mod)) {
            stop("Failed retrieving data for the experiment for the start date ", sdate, ", member ", member, " and chunk ", chunk, ".")
          }
          toto1_exp_member <- abind(toto1_exp_member, toto1_exp_chunk$mod, along = 4)
          if (!('lon' %in% names(toto1))) {
			toto1[['lon']] <- toto1_exp_chunk$lon
			toto1[['lat']] <- toto1_exp_chunk$lat
	      }
        }
        if (load_obs_data) {
          toto1_obs_chunk <- Load(varname, NULL, lstobs, sdate, leadtimemin = ltimes_to_take * (chunk_number - 1) + 1, leadtimemax = ltimes_to_take * chunk_number, output = 'lonlat', grid = grid)
          if (is.null(toto1_obs_chunk$obs)) {
            stop("Failed retrieving data for the observation for the start date ", sdate, ", member ", member, " and chunk ", chunk, ".")
          }
          toto1_obs_sdate <- abind(toto1_obs_sdate, toto1_obs_chunk$obs, along = 4)
        }
      }
	  if (dim(toto1_exp_member)[4] < greatest_n_of_ltimes) {
		padding_array_dims <- dim(toto1_exp_member)
		padding_array_dims[4] <- greatest_n_of_ltimes - padding_array_dims[4]
		toto1_exp_member <- abind(toto1_exp_member, array(dim = padding_array_dims), along = 4)
      }
      toto1_exp_sdate <- abind(toto1_exp_sdate, toto1_exp_member, along = 2)
	  if (dim(toto1_obs_sdate)[4] < greatest_n_of_ltimes) {
		padding_array_dims <- dim(toto1_obs_sdate)
		padding_array_dims[4] <- greatest_n_of_ltimes - padding_array_dims[4]
		toto1_obs_sdate <- abind(toto1_obs_sdate, array(dim = padding_array_dims), along = 4)
      }
	}
	if (!is.null(toto1_exp_sdate)) {
      if (dim(toto1_exp_sdate)[2] < greatest_n_of_members) {
        padding_array_dims <- dim(toto1_exp_sdate)
	    padding_array_dims[4] <- greatest_n_of_ltimes - padding_array_dims[4]
	    toto1_exp_sdate <- abind(toto1_exp_sdate, array(dim = padding_array_dims), along = 4)
      }
      toto1$mod <- abind(toto1$mod, toto1_exp_sdate, along = 3)
	}
	if (!is.null(toto1_obs_sdate)) {
	  toto1$obs <- abind(toto1$obs, toto1_obs_sdate, along = 3)
    }
  }
  if (leading_empty_sdates > 0) {
	padding_array_dims <- dim(toto1$mod)
	padding_array_dims[3] <- leading_empty_sdates
	toto1$mod <- abind(array(dim = leading_empty_sdates), toto1$mod, along = 3)
  }
  if (('b02p' %in% lstexpid)) {
    toto1bis <- Load(varname, 'b02p', obs = NULL, '19501101', output = 'lonlat', grid = grid)
    if (is.null(toto1bis$mod)) {
      stop("Failed retrieving data for b02p.")
    }
    toto1ter <- Histo2Hindcast(toto1bis$mod, '19501101', sdates, 
                               nleadtimesout = ifelse(max(sapply(chunks, length)) == 1, nltimeout, max(sapply(chunks, length)) * nltimechunk))
    toto1beta <- array(dim = c(dim(toto1$mod)[1] + dim(toto1ter)[1], 
                       max(dim(toto1$mod)[2], dim(toto1ter)[2]),
                       dim(toto1$mod)[3:6]))
    toto1beta[1:dim(toto1$mod)[1], 1:dim(toto1$mod)[2], , , , ] <- toto1$mod
    toto1beta[(dim(toto1$mod)[1] + 1):(dim(toto1$mod)[1] + dim(toto1ter)[1]), 
              1:dim(toto1ter)[2], , , , ] <- toto1ter
    toto1$mod <- toto1beta
  }

  toto1$mod <- InsertDim(Mean1Dim(Season(toto1$mod, 4, mon0, mon0, (mon0 + min(11, greatest_n_of_ltimes - 1) - 1) %% 12 + 1), 4), 4, 1)
  toto1$obs <- InsertDim(Mean1Dim(Season(toto1$obs, 4, mon0, mon0, (mon0 + min(11, greatest_n_of_ltimes - 1) - 1) %% 12 + 1), 4), 4, 1)

  if (var == 'pr') {
    toto1$mod <- toto1$mod * 3600 #1000 * 3600 * 24
    toto1$obs <- toto1$obs * 1000 * 3600 * 24
  }
  if (var == 'heatcsum' | var == 'heatcsum1-17' | var == 'heatcsum18-22' | var == 'heatcsum23-46') {
    toto1$mod <- toto1$mod / 1e22
    toto1$obs <- toto1$obs / 1e22
  }
  if (var == 'siarea' | var=='siextent' | var=='sivol') {
    toto1$mod <- toto1$mod/1000
    if (var == 'sivol') {
      toto1$obs <- toto1$obs/1000
    }
  }
  save(toto1, file = paste(savename, '.sav', sep = ''))
}

clims <- Clim(toto1$mod, toto1$obs)
ano_exp <- Ano(toto1$mod, clims$clim_exp)
ano_obs <- Ano(toto1$obs, clims$clim_obs)

if (compcor) {
  cor <- Corr(Mean1Dim(ano_exp, 2), Mean1Dim(ano_obs, 2), 1, 2)
  cols <- c("dodgerblue4", "dodgerblue1", "forestgreen", "yellowgreen",
            "white", "white", "yellow", "orange", "red", "saddlebrown")
  lims <- seq(-1, 1, 0.2)
  for (jexp in 1:length(lstexpid)) {
    flag <- array(F, dim = dim(cor[jexp, 1, 2, 1, , ]))
    flag[which(cor[jexp, 1, 2, 1, , ] > cor[jexp, 1, 4, 1, , ])] <- T
    postscript(paste('CorCoef2d_', var, '_', lstexpid[jexp], '_', season, '_',
               ly, '.eps', sep = ''))
    PlotEquiMap(cor[jexp, 1, 2, 1, , ], toto1$lon, toto1$lat, 
                toptitle = paste('Correlation Coefficient', lstexpid[jexp],
                switch(season, 'Year' = 'annual', season), switch(var, 
                'tas' = 'near surface temperature', 
                'tos' = 'sea surface temperature', 'pr' = 'precipitation'),
                switch(var, 'tas' = 'GHCNv2+ERSSTv3b+GISTEMP', 
                'tas' = 'ERSSTv3b', 'pr' = 'GPCC'), switch(ly, 
                'ly1' = 'Year1', 'ly2-5' = 'Year2-5', 'ly6-9' = 'Year6-9')),
                sizetit = 0.8, brks = lims, cols = cols, colNA = switch(var,
                'pr' = 'white', grey(0.4)), filled.continents = switch(var,
                'tas' = F, 'tos' = T, 'pr' = F), dots = flag, intylat = 45)
    dev.off()
  }
}

if (comprms) {
  rmse <- RMS(Mean1Dim(ano_exp, 2), Mean1Dim(ano_obs, 2), 1, 2)
  cols <- rev(c("dodgerblue4", "dodgerblue1", "forestgreen", "yellowgreen",
                "white", "white", "yellow", "orange", "red", "saddlebrown"))
  lims <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.7, 1, 1.5, 2)
  lims <- switch(var, 'tas' = lims * 2, 'tos' = lims * 2, lims)
  rmse[which(rmse > max(lims))] <- max(lims)
  for (jexp in 1:length(lstexpid)) {
    postscript(paste('RMSE2d_', var, '_', lstexpid[jexp], '_', season, '_', ly,
               '.eps', sep = ''))
    PlotEquiMap(rmse[jexp, 1, 2, 1, , ], toto1$lon, toto1$lat, 
                toptitle = paste('RMSE', lstexpid[jexp], switch(season, 
                'Year' = 'annual', season), switch(var, 
                'tas' = 'near surface temperature', 
                'tos' = 'sea surface temperature', 'pr' = 'precipitation'),
                switch(var, 'tas' = 'GHCNv2+ERSSTv3b+GISTEMP', 
                'tas' = 'ERSSTv3b', 'pr' = 'GPCC'), switch(ly, 
                'ly1' = 'Year1', 'ly2-5' = 'Year2-5', 'ly6-9' = 'Year6-9')),
                sizetit = 0.8, brks = lims, cols = cols, colNA = switch(var, 
                'pr' = 'white', grey(0.4)), filled.continents = switch(var,
                'tas' = F, 'tos' = T, 'pr' = F), intylat = 45)
    dev.off()
  }
}

if (comprmsss) {
  rmsss <- RMSSS(Mean1Dim(ano_exp, 2), Mean1Dim(ano_obs, 2), 1, 2)
  cols <- c("dodgerblue4", "dodgerblue1", "forestgreen", "yellowgreen", 
            "white", "white", "yellow", "orange", "red", "saddlebrown")
  lims <- seq(-1, 1, 0.2)
  for (jexp in 1:length(lstexpid)) {
    flag <- array(F, dim = dim(rmsss[jexp, 1, 2, 1, , ]))
    flag[which(rmsss[jexp, 1, 2, 1, , ] < 0.05)] <- T
    rmsss[which(-1 > rmsss)] = -1
    postscript(paste('RMSSS2d_', var, '_', lstexpid[jexp], '_', season, '_', ly,
               '.eps', sep = ''))
    PlotEquiMap(rmsss[jexp, 1, 1, 1, , ], toto1$lon, toto1$lat, 
                toptitle = paste('RMSSS', lstexpid[jexp], switch(season,
                'Year' = 'annual', season), switch(var, 
                'tas' = 'near surface temperature', 
                'tos' = 'sea surface temperature', 'pr' = 'precipitation'), 
                switch(var, 'tas' = 'GHCNv2+ERSSTv3b+GISTEMP', 
                'tas' = 'ERSSTv3b', 'pr' = 'GPCC'), switch(ly, 
                'ly1' = 'Year1', 'ly2-5' = 'Year2-5', 'ly6-9' = 'Year6-9')),
                sizetit = 0.8, brks = lims, cols = cols, colNA = switch(var,
                'pr' = 'white', grey(0.4)), filled.continents = switch(var,
                'tas' = F, 'tos' = T, 'pr' = F), dots = flag, intylat = 45)
    dev.off()
  }
}

if (compratrms) { 
  ratrms <- RatioRMS(Mean1Dim(ano_exp, 2)[1, , 1, , ], 
                     Mean1Dim(ano_exp, 2)[2, , 1, , ], 
                     Mean1Dim(ano_obs, 2)[1, , 1, , ], 1)
  cols <- rev(c("dodgerblue4", "dodgerblue1", "forestgreen", "yellowgreen",
                "white", "white", "yellow", "orange", "red", "saddlebrown"))
  lims <- c(0, 0.5, 0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2, 2, 6)
  flag <- array(F, dim = dim(ratrms[1, , ]))
  flag[which(ratrms[2, , ] < 0.05)] <- T
  postscript(paste('Rati_RMSE2d_', var, '_', lstexpid[1], '_', lstexpid[2], 
                   '_', season, '_', ly, '.eps', sep = ''))
  PlotEquiMap(ratrms[1, , ], toto1$lon, toto1$lat, toptitle = paste('RMSE',
              lstexpid[1], '/ RMSE', lstexpid[2], switch(season, 
              'Year' = 'annual', season), switch(var, 
              'tas' = 'near surface temperature', 
              'tos' = 'sea surface temperature', 'pr' = 'precipitation'),
              switch(var, 'tas' = 'GHCNv2+ERSSTv3b+GISTEMP', 'tas' = 'ERSSTv3b',
              'pr' = 'GPCC'), switch(ly, 'ly1' = 'Year1', 'ly2-5' = 'Year2-5',
              'ly6-9' = 'Year6-9')), sizetit = 0.8, brks = lims, cols = cols,
              colNA = switch(var, 'pr' = 'white', grey(0.4)), 
              filled.continents = switch(var, 'tas' = F, 'tos' = T, 'pr' = F),
              dots = flag, intylat = 45)
  dev.off()
}
