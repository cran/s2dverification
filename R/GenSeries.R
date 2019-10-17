#'Generates An AR1 Time Series
#'
#'This function generates AR1 processes containing n data points, where alpha 
#'is the autocorrelation at lag 1, and the mean and standard deviation are 
#'specified by the mean and std arguments.
#'
#'@param n Length of the timeseries to be generated.
#'@param alpha Autocorrelation at lag 1.
#'@param mean Mean of the data.
#'@param std Standard deviation of the data.
#'
#'@return AR1 timeseries.
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2012-04  (L. Auger, \email{ludovic.auger@meteo.fr})  -  Original code\cr
#'1.0  -  2012-04  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
#'@examples
#'series <- GenSeries(1000, 0.35, 2, 1)
#'plot(series, type = 'l')
#'
#'@importFrom stats rnorm
#'@export
GenSeries <- function(n, alpha, mean, std) {
  res <- vector("numeric", n)
  x <- mean
  stdterm <- std * (sqrt(1 - alpha ^ 2) / (1 - alpha))
  for (i in 1:100) {
    x <- alpha * x + (1 - alpha) * rnorm(1, mean, stdterm)
  }
  for (i in 1:n) {
    x <- alpha * x + (1 - alpha) * rnorm(1, mean, stdterm)
    res[i] <- x
  }
  
  res
}
