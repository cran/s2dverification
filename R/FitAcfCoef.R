#'Fits an AR1 AutoCorrelation Function Using the Cardano Formula
#'
#'This function finds the minimum point of the fourth order polynom 
#'(a - x)2 + 0.25(b - x2)2 written to fit the two autoregression coefficients 
#'a and b.\cr
#'A consequence of the Cardano formula is that, provided a and b are in [0 1], 
#'the problem is well posed, delta > 0 and there is only one minimum.\cr\cr
#'This function is called in Alpha() to minimize the mean square differences 
#'between the theoretical autocorrelation function of an AR1 and the first 
#'guess of the estimated autocorrelation function estacf, using only the 
#'first two lags.
#' 
#'@param a Coefficient a : first estimate of the autocorrelation at lag 1.
#'@param b Coefficient b : first estimate of the autocorrelation at lag 2.
#'
#'@return Best estimate of the autocorrelation at lag 1.
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2012-06  (L. Auger)  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens)  -  Formatting to CRAN
#'@examples
#'series <- GenSeries(1000, 0.35, 2, 1)
#'estacf <- acf(series[951:1000], plot = FALSE)$acf
#'alpha <- FitAcfCoef(max(estacf[2], 0), max(estacf[3], 0))
#'print(alpha)
#'
#'@export
FitAcfCoef <- function(a, b) {
  if ((a < 0) || (a > 1) || (b < 0) || (b > 1)) {
    stop("One of the coefficients is outside the [0 1] interval");
  }
  p <- (2 - b)
  q <- -2 * a
  delta <- q**2 + 4 / 27 * p**3

  minimum <- (((-q + sqrt(delta)) * 0.5) ^ (1 / 3) 
              - ((q + sqrt(delta)) * 0.5) ^ (1 / 3))
          
  minimum 
}
