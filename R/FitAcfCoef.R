FitAcfCoef <- function(a, b) {
  # This function finds the minimum point of the fourth order polynom 
  # (a - x)2 + 0.25(b - x2)2 written to fit the two autoregression 
  # coefficients a and b.
  # Thanks to the Cardano formula provided a and b in [0 1] the problem is 
  # well posed delta > 0 and there is only one solution to the minimum. 
  #
  # This function is called in Alpha() to minimize the mean square differences 
  # between the theoretical autocorrelation function of an AR1 and the first 
  # guess of estimated autocorrelation function estacf, using only the first 
  # two lags. 
  #
  # Args:
  #   a: Coefficient a.
  #   b: Coefficient b.
  #
  # Returns:
  #   Minimum value of the polynom.
  #
  # History:
  #   1.0  #  2012-06  (L. Auger, ludovic.auger@meteo.fr)  #  Original code

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
