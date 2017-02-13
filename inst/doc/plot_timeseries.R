#!/usr/bin/env Rscript

library(s2dverification)
library(ncdf4)
library(abind)
args <- commandArgs(TRUE)

comptrend <- TRUE
compcor <- TRUE
comprms <- TRUE
compspread <- TRUE
plotano <- TRUE

var <- args[1]    # siextent/siarea/sivol/tos/tas/pr/heatcsum/heatcsum1-17/heatcsum18-22/heatcsum23-46/vsftmyz
pole <- args[2]   # N/S only for siarea/siextent
nltimemax <- 124  # number of leadtimes max in the experiments (in months)
nltimeout <- 60   # number of leadtimes to postprocess(in months)
nltimechunk <- 4  # number of leadtimes per chunk (in months)
lstexpid <- c('i00k', 'b02p') # list of ids
grid <- '320x160' # atmospheric grid for tos/pr (ocean/land only)
vertnem <- 'L42'  # Number of vertical levels in nemo
mon0 <- 11        # initial month
year0 <- 1960     # first start date
yearf <- 2005     # last start date
intsdate <- 5     # interval between start dates
runmeanlen <- 12  # length of the window for running mean (in months)
members <- list(19900101 = c('r1i1p1')) 
#PUT IN ORDER. NONE CAN BE MISSING
chunks <- list(19900101 = c('199001-199001', '199002-199002', '199003-199003'))

obs <- switch(var, 'siarea' = c('hadisst_v1.1'), 'siextent' = c('hadisst_v1.1'), 
              'tas' = c('ncep-reanalysis', 'era40'), 'tos' = c('ersst_v3b', 
                                                               'hadisst_v1.1'),
              'pr' = c('cru_v3.0', 'gpcc1x1_v6'), 'heatcsum' = c('nemovar_system4'),
              'heatcsum18-22' = c('nemovar_system4'), 'heatcsum23-46' = c('nemovar_system4'),
              'heatcsum1-17' = c('nemovar_system4'), 'vsftmyz' = c('nemovar_system4'),
              'sivol' = 'piomas')
toptitle2 <- switch(var, 'siarea' = "sea ice area", 'siextent' = "sea ice extent",
                    'sivol' = "sea ice volume", 'tas' = "global T2m", 
                    'tos' = "global SST (60S-65N)", 
                    'pr' = 'land precipitation (60S-65N)', 
                    'heatcsum' = "global ocean heat content", 
                    'heatcsum1-17' = 'global 800m-bottom ocean heat content',
                    'heatcsum18-22' = 'global 350m-800m ocean heat content',
                    'heatcsum23-46' = 'global 0-350m ocean heat content',
                    'vsftmyz' = 'Atlantic Overturning Streamfunction (40-55N, 1-2km)'
                    )
ytitle1 <- switch(var, 'siarea' = "Millions km2", 'siextent' = "Millions km2", 
                  'sivol' = 'Thousands km3', 'tas' = 'K', 'tos' = 'K', 
                  'pr' = 'mm/day', 'heatcsum' = '10e22 Joules', 
                  'heatcsum1-17' = '10e22 Joules', 'heatcsum18-22' = '10e22 Joules',
                  'heatcsum23-46' = '10e22 Joules', 'vsftmyz' = 'Sv')

syears <- seq(year0, yearf, intsdate)
imon2 <- paste("0", as.character(mon0), sep = "")
sdates <- paste(as.character(syears), substr(imon2, nchar(imon2) - 1, 
                nchar(imon2)), '01', sep = "")
toptitle1 <- paste(switch(pole, 'N' = "Arctic", 'S' = "Antarctic", ""),
                   toptitle2)

savename <- paste(var, switch(pole, 'N' = paste('_', pole, sep = ''),
                  'S' = paste('_', pole, sep = ''), ''), sep = '')
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
  if (var == 'pr' | var == 'tos' ) {
    cat("Retrieving mask files...\n")
    fnc <- nc_open(paste('/esnas/exp/ecearth/constant/land_sea_mask_', grid, '.nc',
                     sep = ''))
    mask <- ncvar_get(fnc, 'LSM')
    nc_close(fnc)
    if (var == 'pr') {
      fnc <- nc_open('/esnas/obs/dwd/gpcc1x1_v6/constant/land_sea_mask.nc'
                       )
      mask_gpcc <- ncvar_get(fnc, 'lsm')
      nc_close(fnc)
      fnc <- nc_open('/esnas/obs/uea/cru_v3.0/constant/mask_cru_land.nc')
      mask_cru <- ncvar_get(fnc, 'pre')
      nc_close(fnc)
      #fnc <- nc_open('/esnas/obs/noaa/gpcp_v2.2/constant_fields/land_sea_mask.nc'
      #                 )
      #mask_gpcp <- ncvar_get(fnc, 'LSM')
      #nc_close(fnc)
      lstmaskobs <- list(mask_cru, mask_gpcc)
    } else {
      mask <- 1 - mask
      lstmaskobs <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
                         NULL)
    }
    lstmask <- list()
    for (iexp in 1:length(lstexpid)) {
       lstmask[[iexp]] <- mask
    }
  } else {
    lstmask <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
                    NULL, NULL, NULL, NULL, NULL) 
    lstmaskobs <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                       NULL)
  }
  latbnd <- switch(var, 'tos' = c(-60, 65), 'pr' = c(-60, 65), c(-90, 90))
  varname <- switch(var, 'siarea' = paste(var, tolower(pole), sep = ''), 
                         'siextent' = paste(var, tolower(pole), sep = ''), 
                         'sivol' = paste(var, tolower(pole), sep = ''), 
                          var)
#  if (is.na(match('b02p', lstexpid)) == TRUE) { 
#    lstload <- lstexpid
#  } else {
#    lstload <- lstexpid[-match('b02p', lstexpid)]
#  }
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
          toto1_exp_chunk <- Load(varname, lstload, NULL, sdate, latmin = latbnd[1],
                                  latmax = latbnd[2], leadtimemax = ltimes_to_take, 
                                  maskmod = lstmask, maskobs = lstmaskobs, dimnames = list(member = 'region'))
          if (is.null(toto1_exp_chunk$mod)) {
            stop("Failed retrieving data for the experiment for the start date ", sdate, ", member ", member, " and chunk ", chunk, ".")
          }
          toto1_exp_member <- abind(toto1_exp_member, toto1_exp_chunk$mod, along = 4)
        }
        if (load_obs_data) {
          toto1_obs_chunk <- Load(varname, NULL, lstobs, sdate, latmin = latbnd[1],
                                  latmax = latbnd[2], 
                                  leadtimemin = ltimes_to_take * (chunk_number - 1) + 1,
                                  leadtimemax = ltimes_to_take * chunk_number, 
                                  maskmod = lstmask, maskobs = lstmaskobs)
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
    toto1bis <- Load(varname, 'b02p', obs = NULL, '19501101', 
                     latmin = latbnd[1], latmax = latbnd[2], maskmod = lstmask,
                     maskobs = lstmaskobs)
    if (is.null(toto1bis$mod)) {
      stop("Failed retrieving data for b02p.")
    }
    toto1ter <- Histo2Hindcast(toto1bis$mod, '19501101', sdates, 
                               nleadtimesout = ifelse(max(sapply(chunks, length)) == 1, nltimeout, max(sapply(chunks, length)) * nltimechunk))
    toto1beta <- array(dim = c(dim(toto1$mod)[1] + dim(toto1ter)[1], 
                       max(dim(toto1$mod)[2], dim(toto1ter)[2]),
                       dim(toto1$mod)[3:4]))
    toto1beta[1:dim(toto1$mod)[1], 1:dim(toto1$mod)[2], , ] <- toto1$mod
    toto1beta[(dim(toto1$mod)[1] + 1):(dim(toto1$mod)[1] + dim(toto1ter)[1]), 
              1:dim(toto1ter)[2], , ] <- toto1ter
    toto1$mod <- toto1beta
  }
  if (var == 'pr') {
    toto1$mod <- toto1$mod * 3600 * 24 #1000 * 3600 * 24
    toto1$obs <- toto1$obs * 1000 * 3600 * 24
  }
  if (var == 'heatcsum' | var == 'heatcsum1-17' | var == 'heatcsum18-22' | var == 'heatcsum23-46') {
    toto1$mod <- toto1$mod / 1e22
    toto1$obs <- toto1$obs / 1e22
  }
  if (var == 'siarea' | var=='siextent' | var=='sivol') {
    #toto1$mod <- toto1$mod/1000
    if (var == 'sivol') {
      toto1$obs <- toto1$obs/1000
    }
  }
  save(toto1, file = paste(savename, '.sav', sep = ''))
}

toto2a <- Clim(toto1$mod, toto1$obs, memb = TRUE)
toto2b_ano_exp <- Ano(toto1$mod, InsertDim(toto2a$clim_exp, 
                                           3, dim(toto1$mod)[3]) )
toto2b_ano_obs <- Ano(toto1$obs, InsertDim(toto2a$clim_obs,
                                           3, dim(toto1$obs)[3]) )
toto3 <- Smoothing(toto2b_ano_exp, runmeanlen, 4)
toto4 <- Smoothing(toto2b_ano_obs, runmeanlen, 4)
suf <- switch(pole, 'N' = paste('_', pole, sep = ''), 'S' = paste('_', pole,
              sep = ''), '')
PlotAno(toto1$mod, toto1$obs, sdates, toptitle = paste(lstexpid, toptitle1),
        ytitle = c(ytitle1, ytitle1, ytitle1), legends = obs, biglab = F, 
        fileout = paste(var, '_', lstexpid, suf, '.eps', sep = ''))
PlotAno(Smoothing(toto1$mod, runmeanlen, 4), 
        Smoothing(toto1$obs, runmeanlen, 4), sdates, 
        toptitle = paste("smoothed", lstexpid, toptitle1),
        ytitle = c(ytitle1, ytitle1, ytitle1), legends = obs, biglab = F, 
        fileout = paste(var, '_', lstexpid, suf, '_smoothed.eps', sep = ''))

if (plotano) {
  PlotAno(toto3, toto4, sdates, toptitle = paste("smoothed", lstexpid,
          toptitle1, "anomalies"), ytitle = c(ytitle1, ytitle1, ytitle1), 
          legends = obs, biglab = F, fileout = paste(var, '_', lstexpid,suf, 
          '_ano.eps', sep = ''))
  PlotClim(toto2a$clim_exp, toto2a$clim_obs, toptitle = paste(switch(pole,
           'N' = "Arctic", 'S' = "Antarctic", ""), toptitle2, "climatologies"),
           ytitle = ytitle1, monini = mon0, listexp = lstexpid, listobs = obs,
           biglab = F, fileout = paste(savename, '_clim.eps', sep = ''))
} 

if (compspread) {
  toto5 <- toto3 - InsertDim(Mean1Dim(toto3, 2, narm = T), 2, dim(toto3)[2])
  toto6 <- Spread(toto5, c(2, 3))
  PlotVsLTime(toto6$iqr, toptitle = paste("InterQuartile Range", toptitle1),
              ytitle = ytitle1, monini = mon0, listexp = lstexpid, biglab = F,
              fileout = paste("IQR_", savename, ".eps", sep = ''))
  PlotVsLTime(toto6$maxmin, toptitle = paste("Maximum-Minimum for", toptitle1),
              ytitle = ytitle1, monini = mon0, listexp = lstexpid, biglab = F, 
              fileout = paste("MaxMin_", savename, ".eps", sep = ''))
  PlotVsLTime(toto6$sd, toptitle = paste("Standard Deviation for", toptitle1),
              ytitle = ytitle1, monini = mon0, listexp = lstexpid, biglab = F,
              fileout = paste("SD_", savename, ".eps", sep = ''))
  PlotVsLTime(toto6$mad, toptitle = paste("Median Absolute Deviation for",
              toptitle1), ytitle = ytitle1, monini = mon0, listexp = lstexpid,
              biglab = F, fileout = paste("Mad_", savename, ".eps", sep = ''))
}

if (compcor) {
  cor <- Corr(Mean1Dim(toto3, 2), Mean1Dim(toto4, 2), 1, 2, compROW = 3, 
              limits = c(ceiling((runmeanlen + 1) / 2), 
              nltimeout - floor(runmeanlen / 2)))
  PlotVsLTime(cor, toptitle = paste("Correlations for", toptitle1), 
              ytitle = "correlation", monini = mon0, limits = c(-1, 2),
              listexp = lstexpid, listobs = obs, biglab = F, 
              hlines = c(-1, 0, 1), fileout = paste("cor_", savename, ".eps",
              sep = ''))
}

if (comprms) {
  rms <- RMS(Mean1Dim(toto3, 2), Mean1Dim(toto4, 2), 1, 2, compROW = 3, 
             limits = c(ceiling((runmeanlen + 1) / 2), 
             nltimeout - floor(runmeanlen / 2)))
  PlotVsLTime(rms, toptitle = paste("RMSE for", toptitle1), ytitle = ytitle1,
              monini = mon0, listexp = lstexpid, listobs = obs, biglab = F,
              fileout = paste("rms_", savename, ".eps", sep = ""))
}

if (comptrend) {
  trends <- Consist_Trend(Mean1Dim(toto3, 2), Mean1Dim(toto4, 2), intsdate / 12)
  PlotVsLTime(trends$trend, toptitle = paste("Trend for", toptitle1), 
              ytitle = paste(ytitle1, "/ year"), monini = mon0,
              listexp = c(lstexpid, obs), biglab = F, fileout = paste("trend_",
              savename, ".eps", sep = ""))
}

rm(list = ls())
quit()
