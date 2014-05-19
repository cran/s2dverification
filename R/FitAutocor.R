FitAutocor <- function(estacf, window = c(-1, 1), prec = 0.01) {
  # This function fits the theoretical autocorrelation function of an AR1 to 
  # the first guess of estimated autocorrelation function estacf containing any
  # number of lags. The fitting relies on a dichotomial minimisation of the 
  # mean square differences between both autocorrelation functions. It returns 
  # the autocorrelation at lag 1 of the fitted AR1 process.
  #
  # Args:
  #   estacf: Estimated autocorrelation function.
  #   window: Epecifies the interval in which lies the autocorrelation at lag 1 
  #           which is sought for.
  #   prec: Determines the precision to which the autocorrelation function at 
  #         lag 1 is to be estimated.
  #
  # Returns:
  #   alpha: Autocorrelation at lag 1 of the fitted AR1 process.
  #
  # History:
  #   1.0  #  2012-02  (V. Guemas, vguemas@ic3.cat)  #  Original code

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
