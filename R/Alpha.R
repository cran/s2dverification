#'Estimates AutoCorrelation At Lag 1 following Guemas et al, BAMS, 2013b
#'
#'This function, relying on the \code{FitAcfCoef()} function, estimates the 
#'autocorrelation at lag 1 of the xdata array following the method described 
#'in Guemas V., Auger L., Doblas-Reyes F., JAMC, 2013. After applying a linear 
#'detrending and/or a filtering of any frequency peak if requested, the sample 
#'autocorrelation is estimated.\cr
#'Then the theoretical autocorrelation of an AR1 is fitted to the sample 
#'autocorrelation using the Cardano's formula (see \code{FitAcfCoef()}) to 
#'obtain the autocorrelation at lag 1. This method assumes xdata is an AR1 
#'process.
#'@param xdata Timeseries from which the autocorrelation at lag 1 is requested.
#'@param detrend TRUE applies a linear detrending to xdata prior to the 
#'  estimation of the autocorrelation at lag 1.
#'@param filter TRUE applies a filtering of any frequency peak prior to the 
#'  estimation of the autocorrelation at lag 1.
#'
#'@return Autocorrelation at lag 1.
#'
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
#'alpha <- Alpha(sampleData$mod[1, 1, , 1])
#'print(alpha)
#'
#'@keywords datagen
#'@author History:\cr
#'  0.1 - 2012-06 (V. Guemas, \email{virginie.guemas@@bsc.es}) - Original code\cr
#'  1.0 - 2013-09 (N. Manubens, \email{nicolau.manubens@@bsc.es}) - Formatting to CRAN
#'@importFrom stats lm confint acf
#'@export
Alpha <- function(xdata, detrend = FALSE, filter = FALSE) {
  tmp <- xdata

  if (detrend == TRUE) {
    reg <- lm(xdata[is.na(xdata) == FALSE] ~ c(1:length(xdata[is.na(xdata) == FALSE])))
    fitted <- array(0, dim = length(xdata[is.na(xdata) == FALSE]))
    if (confint(reg)[2, 1] > 0 & confint(reg)[2, 2] > 0) {
      fitted <- c(1:length(xdata[is.na(xdata) == FALSE])) * min(confint(reg)[2, 2],
                  confint(reg)[2, 1])
    }
    if (confint(reg)[2, 1] < 0 & confint(reg)[2, 2] < 0) {
      fitted <- c(1:length(xdata[is.na(xdata) == FALSE])) * max(confint(reg)[2, 2],
                  confint(reg)[2, 1])
    }
    tmp[is.na(xdata) == FALSE] <- xdata[is.na(xdata) == FALSE] - fitted
  }

  if (filter == TRUE) {
    spec_z <- Spectrum(xdata)
    for (jlen in 1:dim(spec_z)[1]) {
      if (spec_z[jlen, 2] > spec_z[jlen, 4]) {
        tmp <- Filter(tmp, spec_z[jlen, 1])
      }
    }
  }

  estacf <- acf(tmp[is.na(xdata) == FALSE], plot = FALSE)$acf
  ##nacf <- length(estacf)
  ##estacf <- estacf[1:min(3, nacf)]
  ##alpha <- FitAutocor(estacf, c(0, 1), 0.001)
  alpha <- FitAcfCoef(max(estacf[2], 0), max(estacf[3], 0))

  alpha
}
