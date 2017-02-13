RMS <- function(var_exp, var_obs, posloop = 1, posRMS = 2, compROW = NULL, 
                limits = NULL, siglev = 0.95, conf = TRUE) {       
  #
  #  Remove data along compROW dim if there is at least one NA between limits
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (is.null(compROW) == FALSE) {
    if (is.null(limits) == TRUE) {
      limits <- c(1, dim(var_obs)[compROW])
    }
    outrows <- (is.na(Mean1Dim(var_obs, compROW, narm = FALSE, limits))) 
    outrows <- InsertDim(outrows, compROW, dim(var_obs)[compROW])
    var_obs[which(outrows)] <- NA
  }
  
  #
  #  Enlarge var_exp & var_obs to 10 dim + move posloop & posRMS to 1st & 2nd 
  #  pos 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimsvar <- dim(var_exp)
  for (iind in 1:length(dimsvar)) {
    if (iind != posloop & dim(var_obs)[iind] != dimsvar[iind]) { 
      stop("var_exp & var_obs must have same dimensions except along posloop")
    }
  }
  if (dimsvar[posRMS] < 2 ) {
    stop("At least 2 values required to compute RMSE")
  }
  enlvarexp <- Enlarge(var_exp, 10)
  enlvarobs <- Enlarge(var_obs, 10)
  nexp <- dimsvar[posloop]
  nobs <- dim(var_obs)[posloop]
  posaperm <- numeric(10)
  posaperm[1] <- posloop
  posaperm[2] <- posRMS
  posaperm[3:10] <- seq(1, 10)[-c(posloop, posRMS)]
  permvarexp <- aperm(enlvarexp, posaperm)
  permvarobs <- aperm(enlvarobs, posaperm)
  dimsaperm <- dim(permvarexp)
  #
  #  RMS & its confidence interval computation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (conf) {
    nvals <- 3
    dim_rms <- 2
    conf_low <- (1 - siglev) / 2
    conf_high <- 1 - conf_low
  } else {
    nvals <- 1
    dim_rms <- 1
  }
  enlrms <- array(dim = c(nexp, nobs, nvals, dimsaperm[3:10]))
  for (jexp in 1:nexp) {
    for (jobs in 1:nobs) {
      dif <- array(dim = dimsaperm[-1])
      dif[, , , , , , , , ] <- permvarexp[jexp, , , , , , , , 
                                          , ] - permvarobs[jobs, , , , , , , , , ]
      enlrms[jexp, jobs, dim_rms, , , , , , , , ] <- Mean1Dim(dif ** 2, 1, 
                                                              narm = TRUE) ** 0.5
      if (conf) {
        eno <- Eno(dif, 1)
        for (j3 in 1:dimsaperm[3]){
          for (j4 in 1:dimsaperm[4]){
            for (j5 in 1:dimsaperm[5]){
              for (j6 in 1:dimsaperm[6]){
                for (j7 in 1:dimsaperm[7]){
                  for (j8 in 1:dimsaperm[8]){
                    for (j9 in 1:dimsaperm[9]){
                      for (j10 in 1:dimsaperm[10]){
                        ndat <- length(sort(dif[, j3, j4, j5, j6, j7, j8, j9, 
                                                j10]))
                        enlrms[jexp, jobs, 1, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- (eno[j3, j4, j5, j6, j7, j8, j9, 
                                            j10] * enlrms[jexp, jobs, 2, j3, j4, j5, j6, j7, 
                                                          j8, j9, j10] ** 2 / qchisq(conf_high, eno[j3, j4, j5, 
                                                                                                    j6, j7, j8, j9, j10] - 1)) ** 0.5
                        enlrms[jexp, jobs, 3, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- (eno[j3, j4, j5, j6, j7, j8, j9, 
                                            j10] * enlrms[jexp, jobs, 2, j3, j4, j5, j6, j7, 
                                                          j8, j9, j10] ** 2 / qchisq(conf_low, eno[j3, j4, j5, 
                                                                                                   j6, j7, j8, j9, j10] - 1)) ** 0.5
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
  }
  
  dim(enlrms) <- c(nexp, nobs, nvals, dimsvar[-c(posloop, posRMS)])
  #  
  #  Output
  # ~~~~~~~~
  #
  enlrms
}

.RMS <- function(exp, obs, siglev = 0.95, conf = TRUE) {       
  #
  #  RMS & its confidence interval computation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (conf) {
    conf_low <- (1 - siglev) / 2
    conf_high <- 1 - conf_low
  }
  dif <- rowMeans(exp) - obs
  enlrms <- mean(dif ** 2, na.rm = TRUE) ** 0.5
  if (conf) {
    eno <- Eno(dif, 1)
    ndat <- length(sort(dif))
    conf_int <- c((eno * enlrms ** 2 / qchisq(conf_high, eno - 1)) ** 0.5,
                  (eno * enlrms ** 2 / qchisq(conf_low, eno - 1)) ** 0.5)
    names(conf_int) <- c("conf_low","conf_high")
  } else {
    conf_int <- c()
    names(conf_int) <- c() 
  }

  results <- c(enlrms, conf_int) 
  names(results) <- c("rms", names(conf_int))
  return(results)
}
