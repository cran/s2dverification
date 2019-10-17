#'Filter Frequency Peaks From An Array
#'
#'This function filters out the selected frequency from a time series.\cr 
#'The filtering is performed by dichotomy, seeking for a frequency around 
#'the parameter \code{freq} and the phase that maximizes the signal to subtract 
#'from the time series.\cr
#'The maximization of the signal to subtract relies on a minimization of the 
#'mean square differences between the time series (xdata) and the cosine of 
#'the specified frequency and phase.
#'
#'@param xdata Array to be filtered.
#'@param freq Frequency to filter.
#'
#'@return Filtered Array.
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2012-02  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
#'1.0  -  2012-02  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
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
#'@importFrom stats lm
#'@export
Filter <- function(xdata, freq) {
  fac1 <- 1
  fac2 <- 1
  ndat <- length(xdata)
  ndat2 <- length(xdata[is.na(xdata) == FALSE])
  maxi <- 0
  endphase <- 0
  for (jfreq in seq(freq - 0.5 / ndat2, freq + 0.5 / ndat2, 0.1 / (ndat2 * 
       fac1))) {
    for (phase in seq(0, pi, (pi / (10 * fac2)))) {
      xtest <- cos(phase + c(1:ndat) * jfreq * 2 * pi)
      test <- lm(xdata[is.na(xdata) == FALSE] ~ xtest[
              is.na(xdata) == FALSE])$fitted.values
      if (sum(test ^ 2) > maxi) { 
        endphase <- phase
        endfreq <- jfreq
      }
      maxi <- max(sum(test ^ 2), maxi)
    }
  }
  xend <- cos(endphase + c(1:ndat) * endfreq * 2 * pi)
  xdata[is.na(xdata) == FALSE] <- xdata[is.na(xdata) == FALSE] - lm(
                              xdata[is.na(xdata) == FALSE] ~ xend[is.na(xdata) == FALSE]
                              )$fitted.values

  xdata
}
