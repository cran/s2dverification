#'Computes Effective Sample Size Following Guemas et al, BAMS, 2013b
#'
#'This function computes the effective number of independent values in the 
#'xdata array following the method described in 
#'Guemas V., Auger L., Doblas-Reyes F., JAMC, 2013. \code{EnoNew} provides 
#'similar functionality to \code{Eno} but with the added options to remove 
#'the linear trend or filter the frequency.
#'
#'@param xdata A numeric vector.
#'@param detrend Should the linear trend be removed from the data prior to 
#'  the estimation of the equivalent number of independent values.
#'@param filter Should a filtering of the frequency peaks be applied prior 
#'  to the estimation of the equivalent number of independant data.
#'
#'@references
#'Guemas V, Auger L, Doblas-Reyes FJ, Rust H, Ribes A, 2014, Dependencies in 
#'  Statistical Hypothesis Tests for Climate Time Series. Bulletin of the 
#'  American Meteorological Society, 95 (11), 1666-1667.
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2012-06  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
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
#'startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
#'sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'),
#'                                                c('observation'), startDates,
#'                                                leadtimemin = 1,
#'                                                leadtimemax = 4,
#'                                                output = 'lonlat',
#'                                                latmin = 27, latmax = 48,
#'                                                lonmin = -12, lonmax = 40)
#'  }
#'eno <- EnoNew(sampleData$mod[1, 1, , 1, 2, 3])
#'print(eno)
#'
#'@export
EnoNew <- function(xdata, detrend = FALSE, filter = FALSE) {
  alpha <- Alpha(xdata, detrend, filter)
  s <- 0
  n <- length(sort(xdata))
  for (lag in 1:n) {
    s <- s + ((n - lag) / n) * alpha ^ lag
  }

  output <- n / (1 + (2 * s))
}
