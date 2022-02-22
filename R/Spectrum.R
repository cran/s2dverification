#'Estimates Frequency Spectrum
#'
#'This function estimates the frequency spectrum of the xdata array together 
#'with its 95\% and 99\% significance level. The output is provided as an 
#'array with dimensions c(number of frequencies, 4). The column contains the 
#'frequency values, the power, the 95\% significance level and the 99\% one.\cr
#'The spectrum estimation relies on a R built-in function and the significance 
#'levels are estimated by a Monte-Carlo method.
#'
#'@param xdata Array of which the frequency spectrum is required.
#'
#'@return Frequency spectrum with dimensions c(number of frequencies, 4). The 
#'  column contains the frequency values, the power, the 95\% significance 
#'  level and the 99\% one.
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2012-02  (V. Guemas, \email{virginie.guemas@@bsc.es})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@@bsc.es})  -  Formatting to CRAN
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
#'
#'ensmod <- Mean1Dim(sampleData$mod, 2)
#'for (jstartdate in 1:3) {
#'  spectrum <- Spectrum(ensmod[1, jstartdate, ])
#'  for (jlen in 1:dim(spectrum)[1]) {
#'    if (spectrum[jlen, 2] > spectrum[jlen, 4]) {
#'      ensmod[1, jstartdate, ] <- Filter(ensmod[1, jstartdate, ],  
#'                                        spectrum[jlen, 1])
#'    }
#'  }
#'}
#'  \donttest{
#'PlotAno(InsertDim(ensmod, 2, 1), sdates = startDates, fileout =
#'        'filtered_ensemble_mean.eps')
#'  }
#'
#'@importFrom stats spectrum cor rnorm sd quantile
#'@export
Spectrum <- function(xdata) {
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
