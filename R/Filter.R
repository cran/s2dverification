Filter <- function(xdata, freq) {
  # This function filters from the xdata array, the signal of frequency freq. 
  # The filtering is performed by dichotomal seek for the frequency around freq 
  # and the phase that maximizes the signal to subtract to xdata.
  # The maximization of the signal to subtract relies on a minimization of the 
  # mean square differences between xdata and a cosine of given frequency and 
  # phase.  
  # 
  # Params:
  #   xdata: Array of data to be filtered.
  #   freq: Frequency to filter.
  #
  # Returns:
  #   Array of filtered data.
  #
  # History:
  #   1.0  #  2012-02  (V. Guemas, virginie.guemas@ic3.cat)  #  Original code

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
              is.na(xdata) == FALSE])$fitted.value
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
                              )$fitted.value

  xdata
}
