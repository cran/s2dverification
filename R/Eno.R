#'Computes Effective Sample Size With Classical Method
#'
#'Computes the effective number of independent values along the posdim 
#'dimension of a matrix.\cr
#'This effective number of independent observations can be used in 
#'statistical/inference tests.\cr
#'Based on eno function from Caio Coelho from rclim.txt.
#'
#'@param obs Matrix of any number of dimensions up to 10.
#'@param posdim Dimension along which to compute the effective sample size.
#'
#'@return Same dimensions as var but without the posdim dimension.
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-05  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to R CRAN  
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
#'                                                output = 'lonlat', 
#'                                                latmin = 27, latmax = 48, 
#'                                                lonmin = -12, lonmax = 40) 
#'  }
#'sampleData$mod <- Season(sampleData$mod, 4, 11, 1, 12)
#'eno <- Eno(sampleData$mod[1, 1, , 1, , ], 1)
#'PlotEquiMap(eno, sampleData$lon, sampleData$lat)
#'
#'@importFrom stats acf na.pass
#'@export
Eno <- function(obs, posdim) {
  dimsvar <- dim(obs)
  if (is.null(dimsvar)) {
    dimsvar <- length(obs)
  }
  enlobs <- Enlarge(obs, 10)
  outdim <- c(dimsvar, array(1, dim = (10 - length(dimsvar))))
  posaperm <- 1:10
  posaperm[posdim] <- 1
  posaperm[1] <- posdim
  enlobs <- aperm(enlobs, posaperm)
  dimsaperm <- outdim[posaperm]
  #
  #  Loop on all dimensions to compute effective number of observations 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  enleno <- array(dim = c(1, dimsaperm[2:10]))
  for (j2 in 1:dimsaperm[2]) {
    for (j3 in 1:dimsaperm[3]) {
      for (j4 in 1:dimsaperm[4]) {
        for (j5 in 1:dimsaperm[5]) {
          for (j6 in 1:dimsaperm[6]) {
            for (j7 in 1:dimsaperm[7]) {
              for (j8 in 1:dimsaperm[8]) {
                for (j9 in 1:dimsaperm[9]) {
                  for (j10 in 1:dimsaperm[10]) {
                    tmp <- enlobs[, j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    if (length(sort(tmp)) > 1 ) {
                      n <- length(sort(tmp))
                      a <- acf(tmp, lag.max = n - 1, plot = FALSE, 
                               na.action = na.pass)$acf[2:n, 1, 1]
                      s <- 0
                      for (k in 1:(n - 1)) {
                        s <- s + (((n - k) / n) * a[k])
                      }
                      enleno[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- min(n / (1 + (2 * s)), n)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  #
  #  Back to the original dimensions
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #dimsvar <- dimsvar[-posdim]
  if (length(dimsvar) == 1) {
    dimsvar <- 1
  } else {
    dimsvar <- dimsvar[-posdim]
  }
  effnumobs <- array(dim = dimsvar)
  effnumobs[] <- aperm(enleno, posaperm)
  #
  #  Outputs
  # ~~~~~~~~~
  #
  effnumobs
}
