Corr <- function(var_exp, var_obs, posloop = 1, poscor = 2, compROW = NULL, 
                 limits = NULL) {
  # Matrix var_exp & var_obs should have the same dimensions except along 
  # posloop where the length can be different (nexp & nobs).
  # Corr computes correlation for each jexp in 1:nexp and each jobs in 1:nobs 
  # which gives nexp x nobs correlation for each other grid point of the matrix. 
  # The correlations are computed along the poscor dimension. If compROW is 
  # given, the correlations are computed only if rows along the (compROW)th 
  # dimension are complete between limits[1] and limits[2], that  mean with no 
  # NA between limits[1] and limits[2]. 
  # Default: limits[1] = 1 and limits[2] = length(compROW dimension).
  # The confidence interval is computed by a Fisher transformation.
  # The significance level relies on a one-sided student-T distribution.
  #
  # Args:
  #   var_exp: Matrix of experimental data.
  #   var_obs: Matrix of observational data, same dimensions as var_exp except 
  #            along posloop.
  #   posloop: Dimension nobs and nexp.
  #   poscor: Dimension along which correlation are to be computed.
  #   compROW: Data taken into account only if (compROW)th row is complete.
  #            Default = NULL.
  #   limits: Complete between limits[1] & limits[2]. Default = NULL.
  #
  # Returns:
  #   Matrix with c(length(posloop) in var_exp, length(posloop) in var_obs, 4, 
  #   all other dimensions of var_exp & var_obs except poscor).
  #   The third dimension of length 4 corresponds to the lower limit of the 
  #   95% confidence interval, the computed correlation, the upper limit of the 
  #   95% confidence interval and the 95% significance level given by a  
  #   one-sided T-test.
  #
  # History:
  #   1.0  #  2011-04  (V. Guemas, vguemas@ic3.cat)  #  Original code                 
                 
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
  #  Enlarge var_exp & var_obs to 10 dim + move posloop & poscor to 1st & 2nd 
  #  pos 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  dimsvar <- dim(var_exp)
  for (iind in 1:length(dimsvar)) {
    if (iind != posloop & dim(var_obs)[iind] != dimsvar[iind]) { 
      stop("var_exp & var_obs must have same dimensions except along posloop")
    }
  }
  if (dimsvar[poscor] < 3 ) {
    stop("At least 3 values required to compute correlation")
  }
  enlvarexp <- Enlarge(var_exp, 10)
  enlvarobs <- Enlarge(var_obs, 10)
  nexp <- dimsvar[posloop]
  nobs <- dim(var_obs)[posloop]
  posaperm <- numeric(10)
  posaperm[1] <- posloop
  posaperm[2] <- poscor
  posaperm[3:10] <- seq(1, 10)[-c(posloop, poscor)]
  permvarexp <- aperm(enlvarexp, posaperm)
  permvarobs <- aperm(enlvarobs, posaperm)
  dimsaperm <- dim(permvarexp)
  #
  #  Loop to compute correlation for each grid point
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  enlCORR <- array(dim = c(nexp, nobs, 4, dimsaperm[3:10]))
  for (jexp in 1:nexp) {
    for (jobs in 1:nobs) {
      for (j3 in 1:dimsaperm[3]) {
        for (j4 in 1:dimsaperm[4]) {
          for (j5 in 1:dimsaperm[5]) {
            for (j6 in 1:dimsaperm[6]) {
              for (j7 in 1:dimsaperm[7]) {
                for (j8 in 1:dimsaperm[8]) {
                  for (j9 in 1:dimsaperm[9]) {
                    for (j10 in 1:dimsaperm[10]) {
                      tmp1 <- permvarexp[jexp, , j3, j4, j5, j6, j7, j8, j9, 
                                         j10]
                      tmp2 <- permvarobs[jobs, , j3, j4, j5, j6, j7, j8, j9,
                                         j10]
                      if (length(sort(tmp1)) > 0 & length(sort(tmp2)) > 2) {
                        toto <- cor(tmp1, tmp2, use = "pairwise.complete.obs")
                        enlCORR[jexp, jobs, 2, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- toto
                        #eno <- min(Eno(tmp2, 1), Eno(tmp1, 1))
                        eno <- Eno(tmp2, 1)
                        t <- qt(0.95, eno - 2)
                        enlCORR[jexp, jobs, 4, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- sqrt((t * t) / ((t * t) + eno - 2))
                        enlCORR[jexp, jobs, 1, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- tanh(atanh(toto) + qnorm(0.975) / sqrt(
                                        eno - 3))
                        enlCORR[jexp, jobs, 3, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- tanh(atanh(toto) + qnorm(0.025) / sqrt(
                                        eno - 3))
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
  #
  CORR <- array(dim = c(nexp, nobs, 4, dimsvar[-c(posloop, poscor)]))
  CORR[] <- enlCORR
  #  
  #  Output
  # ~~~~~~~~
  #
  CORR
}
