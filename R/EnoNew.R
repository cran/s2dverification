EnoNew <- function(xdata, detrend = FALSE, filter = FALSE) {
  # This function computes the equivalent number of independant data in the 
  # xdata array following the method described in Guemas V., Auger L., 
  # Doblas-Reyes F., JAMC, 2013.
  #
  # Params:
  #   xdata: Array of data.
  #   detrend: T applies a linear detrending to xdata prior to the estimation 
  #            of the number of independant data.
  #   filter: T applies a filtering of any cycle prior to the estimation of the 
  #           number of independant data.
  #
  # Returns:
  #   Equivalent number of independant data.
  #
  # History:
  #   1.0  #  2012-06  (V. Guemas, virginie.guemas@ic3.cat)  #  Original code

  alpha <- Alpha(xdata, detrend, filter)
  s <- 0
  n <- length(sort(xdata))
  for (lag in 1:n) {
    s <- s + ((n - lag) / n) * alpha ^ lag
  }

  output <- n / (1 + (2 * s))
}
