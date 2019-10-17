#'Computes The Regression Of An Array On Another Along A Dimension
#'
#'Computes the regression of the input matrice vary on the input matrice varx 
#'along the posREG dimension by least square fitting. Provides the slope of 
#'the regression, the associated confidence interval, and the intercept.\cr
#'Provides also the vary data filtered out from the regression onto varx.\cr
#'The confidence interval relies on a student-T distribution.
#'
#'@param vary Array of any number of dimensions up to 10.
#'@param varx Array of any number of dimensions up to 10. 
#'  Same dimensions as vary.
#'@param posREG Position along which to compute the regression.
#'
#'@return 
#'  \item{$regression}{
#'    Array with same dimensions as varx and vary except along posREG 
#'    dimension which is replaced by a length 4 dimension, corresponding 
#'    to the lower limit of the 95\% confidence interval, the slope, 
#'    the upper limit of the 95\% confidence interval and the intercept.
#'  }
#'  \item{$filtered}{
#'    Same dimensions as vary filtered out from the regression onto varx 
#'    along the posREG dimension.
#'  }
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2013-05  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
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
#'                   output = 'lonlat', latmin = 27, latmax = 48, 
#'                   lonmin = -12, lonmax = 40)
#'  }
#'  \dontshow{
#'startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
#'sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'),
#'                                                c('observation'), startDates,
#'                                                output = 'lonlat',
#'                                                latmin = 27, latmax = 48,
#'                                                lonmin = -12, lonmax = 40)
#'  }
#'sampleData$mod <- Season(sampleData$mod, 4, 11, 12, 2)
#'sampleData$obs <- Season(sampleData$obs, 4, 11, 12, 2)
#'reg <- Regression(Mean1Dim(sampleData$mod, 2),
#'                  Mean1Dim(sampleData$obs, 2), 2)
#'PlotEquiMap(reg$regression[1, 2, 1, , ], sampleData$lon, sampleData$lat, 
#'            toptitle='Regression of the prediction on the observations', 
#'            sizetit = 0.5)
#'
#'@importFrom stats lm na.omit confint
#'@export
Regression <- function(vary, varx, posREG = 2) {
  #
  #  Enlarge the size of varx and vary to 10
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  dimsvar <- dim(vary)
  if (is.null(dimsvar)) {
    dimsvar <- length(vary)
  }
  if (length(dimsvar) == 1) {
    if (length(varx) != dimsvar[1]) {
      stop("varx and vary should have the same dimensions")
    }
  } else {
    for (jdim in 1:length(dimsvar)) {
      if (dim(varx)[jdim] != dimsvar[jdim]) {
        stop("varx and vary should have the same dimensions")
      }
    }
  }
  enlvarx <- Enlarge(varx, 10)
  enlvary <- Enlarge(vary, 10)
  outdim <- dim(enlvarx)
  
  #
  #  Initialize intermediate and output matrices
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  posaperm <- 1:10
  posaperm[1] <- posREG
  posaperm[posREG] <- 1
  enlvarx <- aperm(enlvarx, posaperm)
  enlvary <- aperm(enlvary, posaperm)
  dimsaperm <- outdim[posaperm]

  enlfilt <- array(dim = dimsaperm)
  filtered <- array(dim = dimsvar)
  dimsaperm[1] <- 4
  enlreg <- array(dim = dimsaperm)
  dimsvar[posREG] <- 4
  regression <- array(dim = dimsvar)
  #
  #  Loop on all dimensions to compute regressions 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  for (j2 in 1:dimsaperm[2]){
    for (j3 in 1:dimsaperm[3]){
      for (j4 in 1:dimsaperm[4]){
        for (j5 in 1:dimsaperm[5]){
          for (j6 in 1:dimsaperm[6]){
            for (j7 in 1:dimsaperm[7]){
              for (j8 in 1:dimsaperm[8]){
                for (j9 in 1:dimsaperm[9]){
                  for (j10 in 1:dimsaperm[10]){
                    tmpy <- enlvary[, j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    tmpx <- enlvarx[, j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    if ((length(sort(tmpy)) > 0) & (length(sort(tmpx)) > 0)) {
                      lm.out <- lm(tmpy ~ tmpx, na.action = na.omit)
                      enlreg[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- confint(lm.out)[2, 1]
                      enlreg[2, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- lm.out$coefficients[2]
                      enlreg[3, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- confint(lm.out)[2, 2]
                      enlreg[4, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- lm.out$coefficients[1]
                      enlfilt[which(is.na(tmpx) == FALSE & is.na(tmpy) == FALSE), j2, 
                              j3, j4, j5, j6, j7, j8, j9, j10] <- tmpy[which(
                              is.na(tmpx) == FALSE & is.na(tmpy
                              ) == FALSE)] - lm.out$fitted.values  
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
  regression[] <- aperm(enlreg, posaperm)
  filtered[] <- aperm(enlfilt, posaperm)

  #
  #  Outputs
  # ~~~~~~~~~
  #
  invisible(list(regression = regression, filtered = filtered))
}
