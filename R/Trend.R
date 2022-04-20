#'Computes the Trend of the Ensemble Mean
#'
#'Computes the trend along the forecast time of the ensemble mean by least 
#'square fitting, and the associated error interval.\cr
#'Trend() also provides the time series of the detrended ensemble mean 
#'forecasts.\cr
#'The confidence interval relies on a student-T distribution.\cr\cr
#'.Trend provides the same functionality but taking a matrix ensemble members 
#'as input (exp).
#'
#'@param var An array of any number of dimensions up to 10.
#'@param posTR An integer indicating the position along which to compute the 
#'  trend.
#'@param interval A number of months/years between 2 points along posTR 
#'  dimension. Set 1 as default.
#'@param siglev A numeric value indicating the confidence level for the 
#'  computation of confidence interval. Set 0.95 as default.
#'@param conf A logical value indicating whether to compute the confidence 
#'  levels or not. Set TRUE as default.
#'@param exp An M by N matrix representing M forecasts from N ensemble members.
#'
#'@return 
#'\item{$trend}{
#'  The intercept and slope coefficients for the least squares fitting of the 
#'  trend. 
#'  An array with same dimensions as parameter 'var' except along the posTR 
#'  dimension, which is replaced by a length 4 (or length 2 if conf = FALSE) 
#'  dimension, corresponding to the lower limit of the confidence interval 
#'  (only present if conf = TRUE), the slope, the upper limit of the confidence 
#'  interval (only present if conf = TRUE), and the intercept.
#'}
#'\item{$detrended}{
#'  Same dimensions as var with linearly detrended var along the posTR 
#'  dimension.
#'}
#'Only in .Trend:
#'\item{$conf.int}{
#'  Corresponding to the limits of the \code{siglev}\% confidence interval 
#'  (only present if \code{conf = TRUE}) for the slope coefficient.
#'} 
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-05  (V. Guemas)  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens)  -  Formatting to CRAN\cr
#'2.0  -  2017-02  (A. Hunter)  -  Adapt to veriApply()
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
#'months_between_startdates <- 60
#'trend <- Trend(sampleData$obs, 3, months_between_startdates)
#'  \donttest{
#'PlotVsLTime(trend$trend, toptitle = "trend", ytitle = "K / (5 year)",
#'            monini = 11, limits = c(-1,1), listexp = c('CMIP5 IC3'),
#'            listobs = c('ERSST'), biglab = FALSE, hlines = 0,
#'            fileout = 'tos_obs_trend.eps')
#'PlotAno(trend$detrended, NULL, startDates, 
#'        toptitle = 'detrended anomalies (along the startdates)', ytitle = 'K',
#'        legends = 'ERSST', biglab = FALSE, fileout = 'tos_detrended_obs.eps')
#'  }
#'
#'@rdname Trend
#'@export
Trend <- function(var, posTR = 2, interval = 1, siglev = 0.95, conf = TRUE) {
  # 
  #  Enlarge the size of var to 10 and move posTR to first position 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
  }
  enlvar <- Enlarge(var, 10)
  outdim <- c(dimsvar, array(1, dim = (10 - length(dimsvar))))
  if (conf) {
    nvals <- 4
    poscoef2 <- 2
    poscoef1 <- 4
  } else {
    nvals <- 2
    poscoef2 <- 1
    poscoef1 <- 2
  }
  outdim[posTR] <- nvals
  posaperm <- 1:10
  posaperm[posTR] <- 1
  posaperm[1] <- posTR
  enlvar <- aperm(enlvar, posaperm)
  dimsaperm <- outdim[posaperm]
  # 
  #  Loop on all dimensions to compute trends 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  enltrend <- array(dim = dimsaperm)
  enldetrend <- array(dim = dim(enlvar))
  for (j2 in 1:dimsaperm[2]) {
    for (j3 in 1:dimsaperm[3]) {
      for (j4 in 1:dimsaperm[4]) {
        for (j5 in 1:dimsaperm[5]) {
          for (j6 in 1:dimsaperm[6]) {
            for (j7 in 1:dimsaperm[7]) {
              for (j8 in 1:dimsaperm[8]) {
                for (j9 in 1:dimsaperm[9]) {
                  for (j10 in 1:dimsaperm[10]) {
                    tmp <- enlvar[, j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    if (any(!is.na(tmp))) {
                      mon <- seq(tmp) * interval
                      lm.out <- lm(tmp ~ mon, na.action = na.omit)
                      enltrend[poscoef2, j2, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- lm.out$coefficients[2]
                      enltrend[poscoef1, j2, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- lm.out$coefficients[1]
                      if (conf) {
                        enltrend[c(1, 3), j2, j3, j4, j5, j6, j7, j8, j9, 
                                 j10] <- confint(lm.out, level = siglev)[2, 1:2]
                      }
                      enldetrend[is.na(tmp) == FALSE, j2, j3, j4, j5, j6, j7, j8, 
                                 j9, j10] <- tmp[is.na(tmp) == FALSE] - lm.out$fitted.values
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
  enldetrend <- aperm(enldetrend, posaperm)
  dim(enldetrend) <- dimsvar 
  enltrend <- aperm(enltrend, posaperm)
  dimsvar[posTR] <- nvals
  dim(enltrend) <- dimsvar
  # 
  #  Outputs
  # ~~~~~~~~~
  # 
  invisible(list(trend = enltrend, detrended = enldetrend))
}

#'@rdname Trend
#'@importFrom stats lm na.omit confint
#'@export
.Trend <- function(exp, interval = 1, siglev = 0.95, conf = TRUE) {
  
  ensmean <- rowMeans(exp, na.rm = TRUE)
  
  if (any(!is.na(ensmean))) {
    mon <- seq(ensmean) * interval
    lm.out <- lm(ensmean ~ mon, na.action = na.omit)
    trend <- c(lm.out$coefficients[2], lm.out$coefficients[1])
    if (conf) {
      conf.int <- confint(lm.out, level = siglev)[2, 1:2]
    }
    detrend <- ensmean[is.na(ensmean) == FALSE] - lm.out$fitted.values
  }  
                 
  # 
  #  Outputs
  # ~~~~~~~~~
  # 
  invisible(list(trend = trend, conf.int = conf.int, detrended = detrend))
}
