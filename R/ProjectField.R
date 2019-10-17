#'Project Anomalies onto Modes of Variability
#'
#'Project anomalies onto modes of variability to get the temporal evolution of 
#'the EOF mode selected. Returns principal components (PCs) by 
#'area-weighted projection onto EOF pattern (from \code{EOF()}). 
#'Able to handle NAs.
#'
#'@param ano Array of forecast or observational reference anomalies from 
#'  \code{Ano()} or \code{Ano_CrossValid} with dimensions (number of forecast 
#'  systems, ensemble members, start dates, forecast horizons, latitudes, 
#'  longitudes).
#'@param eof R object with EOFs from \code{EOF}.
#'@param mode Variability mode number in the provided EOF object which to 
#'  project onto.
#'
#'@return Array of principal components in verification format (number of 
#'  forecast systems, ensemble members, start dates, forecast horizons).
#'
#'@keywords datagen
#'@seealso EOF, NAO, PlotBoxWhisker
#'@author History:\cr
#'0.1  -  2012-03  (F. Lienert, \email{flienert@ic3.cat}  -  Original code\cr
#'0.2  -  2014-03  (Lauriane Batte, \email{lauriane.batte@ic3.cat}  -  Bug-fixes:\cr
#'        1- Extra weighting of the anomalies before projection.\cr
#'        2- Reversion of the anomalies along latitudes.\cr
#'        3- Extra-normalisation not necessary.\cr
#'0.3  -  2014-03  (Virginie Guemas, \email{virginie.guemas@bsc.es}  -  Bug-fixes:\cr
#'        1- Another extra-normalisation.\cr
#'        2- 15 lines to compute the em reduced to 1.
#'0.4  -  2014-03  (Lauriane Batte, \email{lauriane.batte@ic3.cat}  -  Normalization\cr
#'by std before returning PCs to be coherent with EOF().\cr
#'0.5  -  2014-04  (Virginie Guemas, \email{virginie.guemas@bsc.es}  - Fixes:\cr
#'        1- Removal of lon, lat, ncpu and neofs argument unused\cr
#'        2- Security checks ano and eof consistency\cr
#'        3- Removal of the mask which is already contained in the EOFs\cr
#'        4- Removal of the PC normalization since we have chosen in\cr
#'\code{EOF()} to normalize the EOFs and multiply the PCs by the normalization\cr
#'factor and the eigenvalue so that the restitution of the original field is \cr
#'done simply by PC * EOFs\cr
#'        5 - The new convention in \code{EOF()} is to divide by the weights\cr
#'so that the reconstruction of the original field rather than the weighted\cr
#'field is obtained by PC * EOFs. The EOFs need therefore to be multiplied \cr
#'back by the weights before projection so that EOF * t(EOF) = 1\cr
#'        6 - Since W *X = PC * EOF if EOF is multiplied back by the weights,\cr
#'PC = W * X * t(EOF) and X the input field to be projected (X) needs to be\cr
#'multiplied by W. Getting input dimensions.
#'1.0  -  2016-03  (N. Manubens, \email{nicolau.manubens@bsc.es})  -  Formatting to R CRAN\cr
#'                 (J.-P. Baudouin, \email{jean.baudouin@bsc.es})  -  Example code and testing
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
#'                   latmin = 27, latmax = 48, lonmin = -12, lonmax = 40)
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
#'  }
#'# Now ready to compute the EOFs and project.
#'ano <- Ano_CrossValid(sampleData$mod, sampleData$obs)
#'
#'# Compute the EOF of the observation.
#'eof <- EOF(ano$ano_obs[1, 1, , 1, , ], sampleData$lon, sampleData$lat)
#'# check the first mode represent the NAO
#'PlotEquiMap(eof$EOFs[1, , ], sampleData$lon, sampleData$lat, filled.continents = FALSE)
#'
#'mode1_exp <- ProjectField(ano$ano_exp, eof, 1)
#'mode1_obs <- ProjectField(ano$ano_obs, eof, 1)
#'
#'# Plot the forecast and the observation of the first mode
#'# for the last year of forecast
#'plot(mode1_obs[1, 1, dim(sampleData$mod)[3], ], type = "l", ylim = c(-1, 1), lwd = 2)
#'for (i in 1:dim(sampleData$mod)[2]) {
#'  par(new = TRUE)
#'  plot(mode1_exp[1, i, dim(sampleData$mod)[3], ], type = "l", col = rainbow(10)[i], 
#'       ylim = c(-15000, 15000))
#'}
#'
#'@export
ProjectField <- function(ano, eof, mode = 1) {
  # Checking ano
  if (!is.numeric(ano) || !is.array(ano)) {
    stop("Parameter 'ano' must be a numeric array.")
  }
  if (length(dim(ano)) != 6) {
    stop("'ano' must have dimensions c(n. data sets, n. members, n. start dates, n. forecast time steps, n. latitudes, n. longitudes).")
  }

  # Checking eof
  print_error <- FALSE
  if (!is.list(eof)) {
    print_error <- TRUE
  } else {
    if (!all(c('EOFs', 'wght') %in% names(eof))) {
      print_error <- TRUE
    } else {
      if (!is.numeric(eof$EOFs) || !is.array(eof$EOFs)) {
        stop("The component 'EOFs' of parameter 'eof' must be a numeric array.")
      }
      if (length(dim(eof$EOFs)) != 3) {
        stop("The component 'EOFs' of parameter 'eof' must have dimensions c(n. modes, n. latitudes, n. longitudes).")
      }
      if (!is.numeric(eof$wght) || !is.array(eof$wght)) {
        stop("The component 'wght' of parameter 'eof' must be a numeric array.")
      }
      if (length(dim(eof$wght)) != 2) {
        stop("The component 'wght' of parameter 'eof' must have dimensions c(n. latitudes, n. longitudes).")
      }
    }
  }
  if (print_error) {
    stop("Parameter 'eof' must be a list with the components 'EOFs' and 'whgt' as output from EOF().")
  }

  # Checking mode
  if (!is.numeric(mode)) {
    stop("Parameter 'mode' must be numeric.")
  }
  mode <- round(mode)
  if (mode > dim(eof$EOFs)[1]) {
    stop("Parameter 'mode' is greater than the number of available modes in 'eof'.")
  }

  nlon <- dim(ano)[6]
  nlat <- dim(ano)[5]
  nyr <- dim(ano)[3]
  lt <- dim(ano)[4]
  nfc <- dim(ano)[2]
  nmod <- dim(ano)[1]
  
  # Security checks.
  if (dim(eof$EOFs)[2] != nlat || dim(eof$wght)[1] != nlat) {
    stop("Inconsistent number of latitudes between parameter 'eof' and input field")
  }
  if (dim(eof$EOFs)[3] != nlon || dim(eof$wght)[2] != nlon) {
    stop("Inconsistent number of longitudes between parameter 'eof' and input field")
  }
  
  # Initialization of pc.ver.
  pc.ver <- array(NA, dim = c(nmod, nfc, nyr, lt))
  
  # Weigths
  e.1 <- eof$EOFs[mode, , ] * eof$wght
  ano <- ano * InsertDim(InsertDim(InsertDim(InsertDim(eof$wght, 
    1, nmod), 2, nfc), 3, nyr), 4, lt)
  
  for (imo in 1:nmod) {
    for (im in 1:nfc) {
      for (iy in 1:nyr) {
        for (il in 1:lt) {
          if (!all(is.na(ano[imo, im, iy, il, , ]))) {
          pc.ver[imo, im, iy, il] <- sum(ano[imo, im, 
            iy, il, , ] * e.1, na.rm = TRUE)
          }
        }
      }
    }
  }
  
  return(pc.ver)
}

