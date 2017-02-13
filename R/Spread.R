Spread <- function(var, posdim = 2, narm = TRUE, siglev = 0.95, conf = TRUE) {
  # 
  #  Enlarge the size of var to 10 and move all posdim to first position 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  dimsvar <- dim(var)
  enlvar <- Enlarge(var, 10)
  outdim <- c(dimsvar, array(1, dim = (10 - length(dimsvar))))
  posaperm <- 1:10
  posaperm[(length(posdim) + 1):10] <- posaperm[-posdim]
  posaperm[1:length(posdim)] <- posdim
  enlvar <- aperm(enlvar, posaperm)
  dimsaperm <- outdim[posaperm]
  outdim <- 1
  for (jpos in 1:length(posdim)) {
    outdim <- outdim * dimsaperm[jpos]
  }
  outdim <- c(outdim, dimsaperm[(length(posdim) + 1):10])
  enlvar <- array(enlvar, dim = outdim)
  enlvar <- Enlarge(enlvar, 10)
  outdim <- c(outdim, array(1, dim = (10 - length(outdim))))
  if (conf) {
    nvals <- 3
    posvar <- 2
  } else {
    nvals <- 1
    posvar <- 1
  }
  outdim[1] <- nvals
  # 
  #  IQR / Max-Min / STD / Median computation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  iqrvar <- array(dim = outdim)
  maxminvar <- array(dim = outdim)
  sdvar <- array(dim = outdim)
  medvar <- array(dim = outdim)
  if (conf) {
    conf_low <- (1 - siglev) / 2
    conf_high <- 1 - conf_low
  }
  for (j2 in 1:outdim[2]) {
    for (j3 in 1:outdim[3]) {
      for (j4 in 1:outdim[4]) {
        for (j5 in 1:outdim[5]) {
          for (j6 in 1:outdim[6]) {
            for (j7 in 1:outdim[7]) {
              for (j8 in 1:outdim[8]) { 
                for (j9 in 1:outdim[9]) {
                  for (j10 in 1:outdim[10]) {
                    # 
                    #  Compute spread along the posdim dimension 
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    # 
                    whichspread <- enlvar[, j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    iqrvar[posvar, j2, j3, j4, j5, j6, j7, j8, j9, 
                           j10] <- IQR(whichspread, na.rm = narm)
                    maxminvar[posvar, j2, j3, j4, j5, j6, j7, j8, j9, 
                              j10] <- max(whichspread, 
                                      na.rm = narm) - min(whichspread, 
                                                          na.rm = narm)
                    sdvar[posvar, j2, j3, j4, j5, j6, j7, j8, j9, 
                          j10] <- sd(whichspread, na.rm = narm)
                    medvar[posvar, j2, j3, j4, j5, j6, j7, j8, j9, 
                           j10] <- mad(whichspread, na.rm = narm)           
                    if (conf) {
                      # 
                      #  Bootstrapping and recomputing the spread
                      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      # 
                      iqr_bs <- c()
                      maxmin_bs <- c()
                      sd_bs <- c()
                      mad_bs <- c()
                      Num <- length(whichspread)
                      for (jmix in 1:100) {
                        drawings <- round(runif(Num, 0.5, Num + 0.5))
                        iqr_bs <- c(iqr_bs, IQR(whichspread[drawings],  
                                    na.rm = narm))
                        maxmin_bs <- c(maxmin_bs, max(whichspread[drawings],
                                       na.rm=narm) - min(whichspread[drawings],
                                                         na.rm = narm))
                        sd_bs <- c(sd_bs, sd(whichspread[drawings], na.rm = narm))
                        mad_bs <- c(mad_bs, mad(whichspread[drawings], 
                                                na.rm = narm))
                      }
                      #
                      #  95% confidence interval obtained from the bootstrapping 
                      #  results
                      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      #  (This confidence interval is itself quite uncertain but 
                      #   it is not possible to produce miracle with 5 values ...)
                      #
                      iqrvar[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- quantile(iqr_bs, conf_low, na.rm = narm)
                      iqrvar[3, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- quantile(iqr_bs, conf_high, na.rm = narm)
                      maxminvar[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- maxminvar[2, j2, j3, j4, j5, j6, j7, j8, 
                                        j9, j10] + (quantile(maxmin_bs, conf_low, 
                                        na.rm = narm) - quantile(maxmin_bs, 
                                        conf_high, na.rm = narm)) / 2
                      maxminvar[3, j2, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- maxminvar[2, j2, j3, j4, j5, j6, j7, j8, 
                                        j9, j10] - (quantile(maxmin_bs, conf_low, 
                                        na.rm = narm) - quantile(maxmin_bs, 
                                        conf_high, na.rm = narm)) / 2
                      sdvar[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                            j10] <- quantile(sd_bs, conf_low, na.rm = narm)
                      sdvar[3, j2, j3, j4, j5, j6, j7, j8, j9, 
                            j10] <- quantile(sd_bs, conf_high, na.rm = narm)
                      medvar[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- medvar[2, j2, j3, j4, j5, j6, j7, j8, j9, 
                                     j10] + (quantile(mad_bs, conf_low, 
                                     na.rm = narm) - quantile(mad_bs, conf_high, 
                                     na.rm = narm))
                      medvar[3, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- medvar[2, j2, j3, j4, j5, j6, j7, j8, j9, 
                                     j10] - (quantile(mad_bs, conf_low, 
                                     na.rm = narm) - quantile(mad_bs, conf_high, 
                                     na.rm = narm))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  #
  #  Back to the original dimensions
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  posaperm <- posdim[1]
  dimsvar <- dimsvar[-posdim]
  if (posdim[1] > 1) {
    posaperm <- c(posaperm, 1:(posdim[1] - 1)) 
  }
  if (posdim[1] <= length(dimsvar)) {
    posaperm <- c(posaperm, (posdim[1] + 1):(length(posdim) + length(dimsvar)))
  }
  while (length(posaperm) < 10) {
    posaperm <- c(posaperm,length(posaperm) + 1) 
  }
  outdim <- outdim[posaperm]
  outdim <- outdim[1:(length(dimsvar) + 1)]
  iqrvar <- aperm(iqrvar, posaperm)
  dim(iqrvar) <- outdim
  maxminvar <- aperm(maxminvar, posaperm)
  dim(maxminvar) <- outdim
  maxminvar[which(is.infinite(maxminvar))] <- NA
  sdvar <- aperm(sdvar, posaperm)
  dim(sdvar) <- outdim
  medvar <- aperm(medvar, posaperm)
  dim(medvar) <- outdim
  medvar[which(is.infinite(medvar))] <- NA
  # 
  #  Outputs
  # ~~~~~~~~~
  # 
  invisible(list(iqr = iqrvar, maxmin = maxminvar, sd = sdvar, mad = medvar))
}
