#'Computes the North Atlantic Oscillation (NAO) Index
#'
#'Compute the North Atlantic Oscillation (NAO) index based on the leading EOF 
#'of the sea level pressure (SLP) anomalies over the north Atlantic region 
#'(20N-80N, 80W-40E). The PCs are obtained by projecting the forecast and 
#'observed anomalies onto the observed EOF pattern (Pobs) or the forecast 
#'anomalies onto the EOF pattern of the other years of the forecast (Pmod). 
#'By default (ftime_average = 2:4) NAO() computes the NAO index for 1-month 
#'lead seasonal forecasts that can be plotted with BoxPlot(). Returns 
#'cross-validated PCs of the NAO index for forecast (ano_exp) and observations 
#'(ano_obs) based on the leading EOF pattern.
#'
#'@param ano_exp Array of North Atlantic SLP (20N-80N, 80W-40E) forecast 
#'  anomalies from \code{Ano()} or \code{Ano_CrossValid()} with dimensions 
#'  (n. of experimental data sets, n. of ensemble members, n. of start dates, 
#'  n. of forecast time steps, n. of latitudes, n. of longitudes). If only 
#'  NAO of observational data needs to be computed, this parameter can be left 
#'  to NULL (default).
#'@param ano_obs Array of North Atlantic SLP (20N-80N, 80W-40E) observed 
#'  anomalies from \code{Ano()} or \code{Ano_CrossValid()} with dimensions 
#'  (n. of observational data sets, n. of obs. ensemble members, 
#'  n. of start dates, n. of forecast time steps, n. of latitudes, 
#'  n. of longitudes). If only NAO of experimental data needs to be computed, 
#'  this parameter can be left to NULL (default).
#'@param lon Vector with the longitudes of \code{ano_exp} and \code{ano_obs}.
#'@param lat Vector with the latitudes of \code{ano_exp} and \code{ano_obs}.
#'@param ftime_average A vector with the forecast time steps to average across 
#'  defining the target period. Takes by default 2:4, i.e. from 2nd to 4th 
#'  forecast time steps.
#'@param obsproj \code{obsproj = TRUE} will compute the NAO index by 
#'  projecting the forecast anomalies onto the leading EOF of observational 
#'  reference.\cr
#'  \code{obsproj = FALSE} will compute the NAO by first computing the leading 
#'  EOF of the forecast anomalies (in cross-validation mode, i.e. leaving the 
#'  year you are evaluating out), and then projecting forecast anomalies onto 
#'  this EOF.
#'
#'@return 
#'\item{NAO_exp}{
#'  Array of forecast NAO index in verification format (ensemble members, 
#'  start dates).
#'  }
#'\item{NAO_obs}{
#'  Array of observed NAO index in verification format (1, number of start
#'  dates).
#'}
#'\item{EOFs_obs}{
#'  EOFs of the observational references.
#'}
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2013-08  (F. Lienert)  -  Original code\cr
#'0.2  -  2014-03  (V. Guemas)  -  Removing the 
#'  rotation\cr
#'0.3  -  2014-05  (L. Batte)  -  Changes to
#'  simplify function and add Pobs and Pmod options for NAO projection 
#'  calculations\cr
#'0.4  -  2015-03  (L. Batte)  -  Polarity 
#'  check and correction is wrong. Switched to have a negative NAO index when the
#'  anomaly pattern corresponds to NAO-.
#'1.0  -  2016-03  (N. Manubens)  -  
#'  Formatted to CRAN
#'@references
#'Doblas-Reyes, F.J., Pavan, V. and Stephenson, D. (2003). The skill of 
#'  multi-model seasonal forecasts of the wintertime North Atlantic Oscillation. 
#'  Climate Dynamics, 21, 501-514. DOI: 10.1007/s00382-003-0350-4
#'
#'@examples
#'# See examples on Load() to understand the first lines in this example
#'  \dontrun{
#'data_path <- system.file('sample_data', package = 's2dverification')
#'expA <- list(name = 'experiment', path = file.path(data_path,
#'             'model/$EXP_NAME$/$STORE_FREQ$_mean/$VAR_NAME$_3hourly',
#'             '$VAR_NAME$_$START_DATE$.nc'))
#'obsX <- list(name = 'observation', path = file.path(data_path,
#'             '$OBS_NAME$/$STORE_FREQ$_mean/$VAR_NAME$',
#'             '$VAR_NAME$_$YEAR$$MONTH$.nc'))
#'
#'# Now we are ready to use Load().
#'startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
#'sampleData <- Load('tos', list(expA), list(obsX), startDates,
#'                   leadtimemin = 1, leadtimemax = 4, output = 'lonlat',
#'                   latmin = 20, latmax = 90, lonmin = -80, lonmax = 40)
#'  }
#'  \dontshow{
#'startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
#'sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'),
#'                                                c('observation'), startDates,
#'                                                leadtimemin = 1,
#'                                                leadtimemax = 4,
#'                                                output = 'lonlat',
#'                                                latmin = 27, latmax = 48,
#'                                                lonmin = -12, lonmax = 40)
#'# No example data is available over NAO region, so in this example we will 
#'# tweak the available data. In a real use case, one can Load() the data over 
#'# NAO region directly.
#'sampleData$lon[] <- c(40, 280, 340)
#'attr(sampleData$lon, 'first_lon') <- 280
#'attr(sampleData$lon, 'last_lon') <- 40
#'attr(sampleData$lon, 'data_across_gw') <- TRUE
#'sampleData$lat[] <- c(20, 80)
#'attr(sampleData$lat, 'first_lat') <- 20
#'attr(sampleData$lat, 'last_lat') <- 80
#'  }
#'
#'# Now ready to compute the EOFs and project on, for example, the first 
#'# variability mode.
#'ano <- Ano_CrossValid(sampleData$mod, sampleData$obs)
#'# Note that computing the NAO over the region for which there is available 
#'# example data is not the full NAO area: NAO() will raise a warning.
#'nao <- NAO(ano$ano_exp, ano$ano_obs, sampleData$lon, sampleData$lat)
#'# Finally plot the NAO index
#'  \donttest{
#'PlotBoxWhisker(nao$NAO_exp, nao$NAO_obs, "NAO index, DJF", "NAO index (PC1) TOS",
#'        monini = 12, yearini = 1985, freq = 1, "Exp. A", "Obs. X")
#'  }
#'
#'@export
NAO <- function(ano_exp = NULL, ano_obs = NULL, lon, lat, ftime_average = 2:4, obsproj = TRUE) {
  # Checking ano_exp
  if (!is.null(ano_exp)) {
    if (!is.numeric(ano_exp) || !is.array(ano_exp)) {
      stop("Parameter 'ano_exp' must be a numeric array.")
    }
    if (length(dim(ano_exp)) != 6) {
      stop("'ano_exp' must have dimensions c(n. experimental data sets, n. members, n. start dates, n. forecast time steps, n. latitudes, n. longitudes).")
    }
  }

  # Checking ano_obs
  if (!is.null(ano_obs)) {
    if (!is.numeric(ano_obs) || !is.array(ano_obs)) {
      stop("Parameter 'ano_obs' must be a numeric array.")
    }
    if (length(dim(ano_obs)) != 6) {
      stop("'ano_obs' must have dimensions c(n. observational data sets, n. obs. members, n. start dates, n. forecast time steps, n. latitudes, n. longitudes).")
    }
  }

  # Checking consistency in ano_exp and ano_obs
  if (!is.null(ano_obs) && !is.null(ano_exp)) {
    if (!identical(dim(ano_exp)[3:6], dim(ano_obs)[3:6])) {
      stop("'ano_obs' and 'ano_exp' must have the same number of start dates, forecast time steps, latitudes and longitudes.")
    }
  }

  # Checking lon and lat
  if (!is.numeric(lon) || !is.numeric(lat)) {
    stop("'lon' and 'lat' must be numeric vectors.")
  }
  if (is.null(attr(lon, 'first_lon')) || is.null(attr(lon, 'last_lon')) || 
      is.null(attr(lon, 'array_across_gw')) || is.null(attr(lon, 'data_across_gw'))) {
    .warning("At least one of the attributes 'first_lon', 'last_lon', 'data_across_gw' or 'array_across_gw' of the parameter 'lon' is not defined (see documentation on output 'lon' of ?Load). The spatial domain of the provided data may be unnoticedly wrong.")
  }
  if (is.null(attr(lat, 'last_lat'))) {
    attr(lat, 'last_lat') <- tail(lat, 1)
  }
  if (is.null(attr(lat, 'first_lat'))) {
    attr(lat, 'first_lat') <- head(lat, 1)
  }
  stop_bad_domain <- "The typical domain used to compute the NAO is 20N-80N, 80W-40E.\n"
  stop_needed <- FALSE
  if (attr(lat, 'last_lat') < 70 || attr(lat, 'last_lat') > 90 || 
      attr(lat, 'first_lat') > 30 || attr(lat, 'first_lat') < 10) {
    stop_needed <- TRUE
  }
  if (!is.null(attr(lon, 'data_across_gw'))) {
    if (!attr(lon, 'data_across_gw')) {
      stop_needed <- TRUE
    }
  }
  if (!is.null(attr(lon, 'first_lon')) && !is.null(attr(lon, 'last_lon'))) {
    if (!(attr(lon, 'last_lon') < attr(lon, 'first_lon'))) {
      stop_needed <- TRUE
    } else if (attr(lon, 'last_lon') > 50 || attr(lon, 'last_lon') < 30 ||
               attr(lon, 'first_lon') > 290 || attr(lon, 'first_lon') < 270) {
      stop_needed <- TRUE
    }
  }
  if (stop_needed) {
    stop(stop_bad_domain)
  }
  ##if (!is.null(attr(lon, 'array_across_gw'))) {
  ##  if (!attr(lon, 'array_across_gw')) {
  ##    REORDER
  ##  }
  ##}

  if (!is.null(ano_exp)) {
    dims <- dim(ano_exp)
  } else if (!is.null(ano_obs)) {
    dims <- dim(ano_obs)
  } else {
    stop("Either one of 'ano_exp' or 'ano_obs' must be provided.")
  }
  nlon <- dims[6]
  nlat <- dims[5]
  nftimes <- dims[4]
  nyr <- dims[3]
  nmemb <- dims[2]
  nexp <- dims[1]

  # Security check lon & lat
  if (length(lon) != nlon) {
    stop("Inconsistent number of longitudes and input field dimensions.")
  }
  if (length(lat) != nlat) {
    stop("Inconsistent number of latitudes and input field dimensions.")
  }

  # Security check nyr
  if (nyr < 2) {
    stop("At least data for 2 start dates must be provided.")
  }

  # Checking ftime_average
  if (!is.numeric(ftime_average)) {
    stop("'ftime_average' must be a numeric vector.")
  }
  if (any(ftime_average > nftimes)) {
    stop("'ftime_averages' contains indexes to non-existing forecast time steps.")
  }

  # Checking obsproj
  if (!is.logical(obsproj)) {
    stop("'obsproj' must be either TRUE or FALSE.")
  }
  if (obsproj) {
    if (is.null(ano_obs)) {
      stop("Parameter 'obsproj' set to TRUE but no 'ano_obs' provided.")
    }
    if (is.null(ano_exp)) {
      .warning("parameter 'obsproj' set to TRUE but no 'ano_exp' provided.")
    }
  }

  fcsys <- 1

  NAOF.ver <- NULL
  NAOO.ver <- NULL
  OEOF <- NULL
  if (!is.null(ano_exp)) {
    ## Target period mean
    ano_exp <- ano_exp[fcsys, , , ftime_average, , , drop = FALSE]
    f1 <- Mean1Dim(ano_exp, posdim = 4, narm = TRUE)
    dim(f1) <- c(1, nmemb, nyr, 1, nlat, nlon)
    ## Cross-validated PCs. Fabian. This should be extended to
    ## nmod and nlt by simple loops. Virginie
    NAOF.ver <- array(NA, c(nmemb, nyr))
  }

  if (!is.null(ano_obs)) {
    ano_obs <- ano_obs[1, 1, , ftime_average, , , drop = FALSE]
    ##gc()
    o1 <- Mean1Dim(ano_obs, posdim = 4, narm = TRUE)
    dim(o1) <- c(1, 1, nyr, 1, nlat, nlon)
    NAOO.ver <- array(NA, c(1, nyr))
  }

  for (iy in 1:nyr) {
    if (!is.null(ano_obs)) {
      ## Observed EOF excluding one forecast start year.
      o2 <- o1[1, 1, -iy, 1, , ]
      dim(o2) <- c(nyr - 1, nlat, nlon)
      OEOF <- EOF(o2, lon, lat, neofs = 1)
      ## Correct polarity of pattern.
      sign <- 1     
      if (0 < mean(OEOF$EOFs[1, which.min(abs(lat - 65)), ], na.rm = T)) {
        sign <- -1
      }
      OEOF$EOFs <- OEOF$EOFs * sign
      OEOF$PCs <- OEOF$PCs * sign
      ## Project observed anomalies.
      PCO <- ProjectField(o1, OEOF, mode = 1)
      ## Keep PCs of excluded forecast start year. Fabian.
      NAOO.ver[1, iy] <- PCO[1, 1, iy, 1]
    }

    if (!is.null(ano_exp)) {
      if (!obsproj) {
        f2 <- f1[1, , -iy, 1, , ]
        dim(f2) <- c(nmemb * (nyr - 1), nlat, nlon)
        FEOF <- EOF(f2, lon, lat, neofs = 1)
        ## Correct polarity of pattern.
        sign <- 1
        if (0 < FEOF$EOFs[1, which.min(abs(lat - 65)), ]) {
          sign <- -1
        }
        FEOF$EOFs <- FEOF$EOFs * sign
        FEOF$PCs <- FEOF$PCs * sign
        ### Lines below could be simplified further by computing
        ### ProjectField() only on the year of interest... (though this is
        ### not vital). Lauriane
        PCF <- ProjectField(f1, FEOF, mode = 1)
        for (imemb in 1:nmemb) {
          NAOF.ver[imemb, iy] <- PCF[1, imemb, iy, 1]
        }
      } else {
        ## Project forecast anomalies on obs EOF
        PCF <- ProjectField(f1, OEOF, mode = 1)
        NAOF.ver[, iy] <- PCF[1, , iy, 1]
      }
    }
  }

  return(list(NAO_exp = NAOF.ver, NAO_obs = NAOO.ver, EOFs_obs = OEOF))
}
