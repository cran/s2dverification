#'Fits an AR1 Autocorrelation Function Using Dichotomy
#'
#'This function fits the theoretical autocorrelation function of an AR1 to 
#'the first guess of the estimated autocorrelation function estacf containing 
#'any number of lags. The fitting relies on a dichotomial minimisation of the 
#'mean square differences between both autocorrelation functions. It returns 
#'the autocorrelation at lag 1 of the fitted AR1 process.
#'
#'@param estacf First guess for the autocorrelation function.
#'@param window Interval in which the autocorrelation at lag 1 should be found.
#'@param prec Precision to which the autocorrelation function at lag 1 is to be estimated.
#'
#'@return Best estimate of the autocorrelation at lag 1.
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2012-02  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
#'@examples
#'series <- GenSeries(1000, 0.35, 2, 1)
#'estacf <- acf(series[951:1000], plot = FALSE)$acf
#'alpha <-  FitAutocor(estacf, c(-1, 1), 0.01)
#'print(alpha)
#'@export
FitAutocor <- function(estacf, window = c(-1, 1), prec = 0.01) {
  nacf <- length(estacf)
  alpha <- mean(window)
  dist <- nacf * 4
  for (jind in 1:nacf) {
    estacf[jind] <- max(estacf[jind], 0)
  }
  for (jalph in seq(window[1], window[2], prec)) {
    test <- sum(((estacf - (jalph ^ c(1:nacf - 1))) / c(1:nacf)) ^ 2)
    if (test < dist) {
      dist <- test
      alpha <- jalph
    }
  }
  alpha
}
