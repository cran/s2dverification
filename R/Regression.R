Regression <- function(vary, varx, posREG = 2) {
  # Compute the regression of matrice vary on matrice varx along the (posREG)th 
  # dimension by least square fitting. Provides the slope of the regression, 
  # the associated confidence interval, and the intercept.
  # Provide also the vary data filtered out from the regression onto varx.
  # The confidence interval relies on a student-T distribution.
  #
  # Args:
  #   varx: Matrix of any number of dimensions up to 10.
  #   vary: Matrix of any number of dimensions up to 10. 
  #   posREG: Position along which to compute the regression.
  #
  # Returns:
  #   $regression: Matrix with same dimensions as varx and vary except along 
  #                posREG dimension which is replaced by a length 4 dimension, 
  #                corresponding to the lower limit of the  95% confidence 
  #                interval, the computed slope and the upper limit of the 95% 
  #                confidence interval and the intercept for each point of the 
  #                matrice along all the other dimensions.
  #   $filtered: Same dimensions as vary filtered out from the regression onto 
  #              varx along the posREG dimension.
  #
  # History:
  #   1.0  #  2013-05  (V. Guemas, vguemas@ic3.cat)  #  Original code   

  #
  #  Enlarge the size of varx and vary to 10
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  dimsvar <- dim(vary)
  if (is.null(dimsvar)) {
    dimsvar <- length(vary)
  }
  if (length(dimsvar) == 1) {
    if (length(varx) != dimsvar[1]) {
      stop("varx and vary should have the same dimensions")
    }
  } else {
    for (jdim in 1:length(dimsvar)) {
      if (dim(varx)[jdim] != dimsvar[jdim]) {
        stop("varx and vary should have the same dimensions")
      }
    }
  }
  enlvarx <- Enlarge(varx, 10)
  enlvary <- Enlarge(vary, 10)
  outdim <- dim(enlvarx)
  
  #
  #  Initialize intermediate and output matrices
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  posaperm <- 1:10
  posaperm[1] <- posREG
  posaperm[posREG] <- 1
  enlvarx <- aperm(enlvarx, posaperm)
  enlvary <- aperm(enlvary, posaperm)
  dimsaperm <- outdim[posaperm]

  enlfilt <- array(dim = dimsaperm)
  filtered <- array(dim = dimsvar)
  dimsaperm[1] <- 4
  enlreg <- array(dim = dimsaperm)
  dimsvar[posREG] <- 4
  regression <- array(dim = dimsvar)
  #
  #  Loop on all dimensions to compute regressions 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  for (j2 in 1:dimsaperm[2]){
    for (j3 in 1:dimsaperm[3]){
      for (j4 in 1:dimsaperm[4]){
        for (j5 in 1:dimsaperm[5]){
          for (j6 in 1:dimsaperm[6]){
            for (j7 in 1:dimsaperm[7]){
              for (j8 in 1:dimsaperm[8]){
                for (j9 in 1:dimsaperm[9]){
                  for (j10 in 1:dimsaperm[10]){
                    tmpy <- enlvary[, j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    tmpx <- enlvarx[, j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    if ((length(sort(tmpy)) > 0) & (length(sort(tmpx)) > 0)) {
                      lm.out <- lm(tmpy ~ tmpx, na.action = na.omit)
                      enlreg[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- confint(lm.out)[2, 1]
                      enlreg[2, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- lm.out$coefficients[2]
                      enlreg[3, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- confint(lm.out)[2, 2]
                      enlreg[4, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- lm.out$coefficients[1]
                      enlfilt[which(is.na(tmpx) == FALSE & is.na(tmpy) == FALSE), j2, 
                              j3, j4, j5, j6, j7, j8, j9, j10] <- tmpy[which(
                              is.na(tmpx) == FALSE & is.na(tmpy
                              ) == FALSE)] - lm.out$fitted.values  
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
  regression[] <- aperm(enlreg, posaperm)
  filtered[] <- aperm(enlfilt, posaperm)

  #
  #  Outputs
  # ~~~~~~~~~
  #
  invisible(list(regression = regression, filtered = filtered))
}
