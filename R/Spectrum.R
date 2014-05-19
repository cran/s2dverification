Spectrum <- function(xdata) {
  # This function estimates the frequency spectrum of the xdata array together 
  # with its 95% and 99% significance level. The output is provided as a matrix 
  # with dimensions c(number of frequencies, 4). The column contains the
  # frequency values, the power, the 95% significance level and the 99% one. 
  # The spectrum estimation relies on a R built-in function and the 
  # significance levels are estimated by a Monte-Carlo method.
  #
  # Params:
  #   xdata: Input array of data.
  #
  # Returns:
  #   The frequency spectrum.
  #
  # History:
  #   1.0  #  2012-02  (V. Guemas, virginie.guemas@ic3.cat)  #  Original code

  print('Warning : Your data are assumed to be evenly spaced in time')
  xdata <- xdata[is.na(xdata) == FALSE]
  ndat <- length(xdata)

  if (ndat >= 3) {
    tmp <- spectrum(xdata, plot = FALSE)
    output <- array(dim = c(length(tmp$spec), 4))
    output[, 1] <- tmp$freq
    output[, 2] <- tmp$spec
    ntir <- 100
    store <- array(dim = c(ntir, length(tmp$spec)))
    for (jt in 1:ntir) {
      toto <- mean(xdata)
      alpha1 <- cor(xdata[2:ndat], xdata[1:(ndat - 1)])
      for (ind in 2:ndat) { 
        b <- rnorm(1, mean(xdata) * (1 - alpha1), sd(xdata) * sqrt(1 - 
                   alpha1 ^ 2))
        toto <- c(toto, toto[ind - 1] * alpha1 + b)
      }
      toto2 <- spectrum(toto, plot = FALSE)
      store[jt, ] <- toto2$spec
    }
    for (jx in 1:length(tmp$spec)) {
      output[jx, 3] <- quantile(store[, jx], 0.95)
      output[jx, 4] <- quantile(store[, jx], 0.99)
    }
  } else {
    output <- NA
  }
  
  output
}
