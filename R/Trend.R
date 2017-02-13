Trend <- function(var, posTR = 2, interval = 1, siglev = 0.95, conf = TRUE) {
  # 
  #  Enlarge the size of var to 10 and move posTR to first position 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
  }
  enlvar <- Enlarge(var, 10)
  outdim <- c(dimsvar, array(1, dim = (10 - length(dimsvar))))
  if (conf) {
    nvals <- 4
    poscoef2 <- 2
    poscoef1 <- 4
  } else {
    nvals <- 2
    poscoef2 <- 1
    poscoef1 <- 2
  }
  outdim[posTR] <- nvals
  posaperm <- 1:10
  posaperm[posTR] <- 1
  posaperm[1] <- posTR
  enlvar <- aperm(enlvar, posaperm)
  dimsaperm <- outdim[posaperm]
  # 
  #  Loop on all dimensions to compute trends 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  enltrend <- array(dim = dimsaperm)
  enldetrend <- array(dim = dim(enlvar))
  for (j2 in 1:dimsaperm[2]) {
    for (j3 in 1:dimsaperm[3]) {
      for (j4 in 1:dimsaperm[4]) {
        for (j5 in 1:dimsaperm[5]) {
          for (j6 in 1:dimsaperm[6]) {
            for (j7 in 1:dimsaperm[7]) {
              for (j8 in 1:dimsaperm[8]) {
                for (j9 in 1:dimsaperm[9]) {
                  for (j10 in 1:dimsaperm[10]) {
                    tmp <- enlvar[, j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    if (any(!is.na(tmp))) {
                      mon <- seq(tmp) * interval
                      lm.out <- lm(tmp ~ mon, na.action = na.omit)
                      enltrend[poscoef2, j2, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- lm.out$coefficients[2]
                      enltrend[poscoef1, j2, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- lm.out$coefficients[1]
                      if (conf) {
                        enltrend[c(1, 3), j2, j3, j4, j5, j6, j7, j8, j9, 
                                 j10] <- confint(lm.out, level = siglev)[2, 1:2]
                      }
                      enldetrend[is.na(tmp) == FALSE, j2, j3, j4, j5, j6, j7, j8, 
                                 j9, j10] <- tmp[is.na(tmp) == FALSE] - lm.out$fitted.values
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
  enldetrend <- aperm(enldetrend, posaperm)
  dim(enldetrend) <- dimsvar 
  enltrend <- aperm(enltrend, posaperm)
  dimsvar[posTR] <- nvals
  dim(enltrend) <- dimsvar
  # 
  #  Outputs
  # ~~~~~~~~~
  # 
  invisible(list(trend = enltrend, detrended = enldetrend))
}

.Trend <- function(exp, interval = 1, siglev = 0.95, conf = TRUE) {
  
  ensmean <- rowMeans(exp, na.rm = TRUE)
  
  if (any(!is.na(ensmean))) {
    mon <- seq(ensmean) * interval
    lm.out <- lm(ensmean ~ mon, na.action = na.omit)
    trend <- c(lm.out$coefficients[2], lm.out$coefficients[1])
    if (conf) {
      conf.int <- confint(lm.out, level = siglev)[2, 1:2]
    }
    detrend <- ensmean[is.na(ensmean) == FALSE] - lm.out$fitted.values
  }  
                 
  # 
  #  Outputs
  # ~~~~~~~~~
  # 
  invisible(list(trend = trend, conf.int = conf.int, detrended = detrend))
}
