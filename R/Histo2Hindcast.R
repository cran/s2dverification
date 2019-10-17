#'Chunks Long Simulations For Comparison With Hindcasts
#'
#'This function reorganizes a long run (historical typically) with only one 
#'start date into chunks corresponding to a set of start dates. The expected 
#'input structure is the one output from \code{Load()} with 4 to 7 dimensions.
#'
#'@param varin Array of model or observational data with dimensions:\cr
#'  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltimes) up to\cr
#'  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltimes, nlevel, nlat, nlon)
#'@param sdatesin Start date of the input matrix 'YYYYMMDD'.
#'@param sdatesout List of start dates of the output matrix 
#'  c('YYYYMMDD', 'YYYYMMDD', ...).
#'@param nleadtimesout Number of leadtimes in the output matrix.
#'
#'@return An array with the same number of dimensions as varin, the same 
#' dimensions 1 and 2 and potentially the same dimensions 5 to 7. Dimensions 
#' 3 and 4 are set by the arguments sdatesout and nleadtimesout.
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2012-11  (V. Guemas, \email{vguemas@ic3.cat})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
#'@examples
#'# See examples on Load() to understand the first lines in this example
#'  \dontrun{
#'data_path <- system.file('sample_data', package = 's2dverification')
#'exp <- list(
#'         name = 'experiment',
#'         path = file.path(data_path, 'model/$EXP_NAME$/monthly_mean',
#'                          '$VAR_NAME$_3hourly/$VAR_NAME$_$START_DATES$.nc')
#'       )
#'obs <- list(
#'         name = 'observation',
#'         path = file.path(data_path, 'observation/$OBS_NAME$/monthly_mean',
#'                          '$VAR_NAME$/$VAR_NAME$_$YEAR$$MONTH$.nc')
#'       )
#'# Now we are ready to use Load().
#'startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
#'sampleData <- Load('tos', list(exp), list(obs), startDates,
#'                   leadtimemin = 1, leadtimemax = 4, output = 'lonlat',
#'                   latmin = 27, latmax = 48, lonmin = -12, lonmax = 40)
#'  }
#'  \dontshow{
#'startDates <- c('19901101')
#'sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'),
#'                                                c('observation'), startDates,
#'                                                leadtimemin = 1,
#'                                                leadtimemax = 4,
#'                                                output = 'areave',
#'                                                latmin = 27, latmax = 48,
#'                                                lonmin = -12, lonmax = 40)
#'  }
#'
#'
#'start_dates_out <- c('19901101', '19911101', '19921101', '19931101', '19941101')
#'leadtimes_per_startdate <- 12
#'experimental_data <- Histo2Hindcast(sampleData$mod, startDates[1], 
#'                                    start_dates_out, leadtimes_per_startdate)
#'observational_data <- Histo2Hindcast(sampleData$obs, startDates[1], 
#'                                     start_dates_out, leadtimes_per_startdate)
#'  \donttest{
#'PlotAno(experimental_data, observational_data, start_dates_out, 
#'        toptitle = paste('anomalies reorganized into shorter chunks'), 
#'        ytitle = 'K', fileout='tos_histo2hindcast.eps')
#'  }
#'
#'@export
Histo2Hindcast <- function(varin, sdatesin, sdatesout, nleadtimesout) {
  #
  #  Input parameters
  # ~~~~~~~~~~~~~~~~~~
  #
  ndims <- length(dim(varin))
  varin <- Enlarge(varin, 7)
  outvar <- array(dim = c(dim(varin)[1:2], length(sdatesout), nleadtimesout,
                  dim(varin)[5:7]))
  yearsin <- as.numeric(substr(sdatesin, 1, 4))
  yearsout <- as.numeric(substr(sdatesout, 1, 4))
  monthin <- as.numeric(substr(sdatesin, 5, 6))
  monthout <- as.numeric(substr(sdatesout, 5, 6))
  #
  #  Re-organization
  # ~~~~~~~~~~~~~~~~~
  #
  for (indyear in 1:length(sdatesout)) {
    difmonths <- (yearsin[1] - yearsout[indyear]
                 ) * 12 + monthin[1] - monthout[1]
    if ((difmonths < nleadtimesout) & ((-difmonths) < dim(varin)[4])) {
      outvar[, , indyear, max(difmonths + 1, 1):min(dim(varin)[4] + difmonths,
             nleadtimesout), , , ] <- varin[, , 1, 
                                        max(1 - difmonths, 
                                        1):min(nleadtimesout - difmonths, 
                                        dim(varin)[4]), , , ]
    }
  }
  #
  #  Outputs
  # ~~~~~~~~~
  #
  outvar <- array(outvar, dim = dim(outvar)[1:ndims])
}
