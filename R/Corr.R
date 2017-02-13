Corr <- function(var_exp, var_obs, posloop = 1, poscor = 2, compROW = NULL, 
                 limits = NULL, siglev = 0.95, method = 'pearson', 
                 conf = TRUE, pval = TRUE) {
  #
  #  Remove data along compROW dim if there is at least one NA between limits
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (!is.null(compROW)) {
    if (is.null(limits)) {
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
  if (method != "kendall" && method != "spearman" && method != "pearson") {
    stop("Wrong correlation method")
  }
  nexp <- dimsvar[posloop]
  nobs <- dim(var_obs)[posloop]
  var_exp <- Enlarge(var_exp, 10)
  var_obs <- Enlarge(var_obs, 10)
  posaperm <- numeric(10)
  posaperm[1] <- posloop
  posaperm[2] <- poscor
  posaperm[3:10] <- seq(1, 10)[-c(posloop, poscor)]
  var_exp <- aperm(var_exp, posaperm)
  var_obs <- aperm(var_obs, posaperm)
  dimsaperm <- dim(var_exp)
  #
  
  # Check the siglev arguments:
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (siglev > 1 || siglev < 0) {
    stop("siglev need to be higher than O and lower than 1")
  }
  #
  #  Loop to compute correlation for each grid point
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dim_val <- 2
  dim_pval <- 4
  nvals <- 1 + 2*conf + pval
  if (!conf) {
    dim_val <- 1
    dim_pval <- 2
  } else {
    conf_low <- (1 - siglev) / 2
    conf_high <- 1 - conf_low
  }
  CORR <- array(dim = c(nexp, nobs, nvals, dimsaperm[3:10]))
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
                      tmp1 <- var_exp[jexp, , j3, j4, j5, j6, j7, j8, j9, 
                                      j10]
                      tmp2 <- var_obs[jobs, , j3, j4, j5, j6, j7, j8, j9,
                                      j10]
                      if (any(!is.na(tmp1)) && sum(!is.na(tmp2)) > 2) {
                        toto <- cor(tmp1, tmp2, use = "pairwise.complete.obs", method = method)
                        CORR[jexp, jobs, dim_val, j3, j4, j5, j6, j7, j8, j9, j10] <- toto
                        #eno <- min(Eno(tmp2, 1), Eno(tmp1, 1))
                        if (pval || conf) {
                          if (method == "kendall" | method == "spearman") {
                            eno <- Eno(rank(tmp2), 1)                          
                          } else if (method == "pearson") {
                            eno <- Eno(tmp2, 1)                          
                          }
                        }
                        if (pval) {
                          #t <- qt(0.95, eno - 2)
                          t <- qt(siglev, eno - 2)
                          CORR[jexp, jobs, dim_pval, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- sqrt((t * t) / ((t * t) + eno - 2))
                        }
                        if (conf) {
                          CORR[jexp, jobs, 1, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- tanh(atanh(toto) + qnorm(conf_high) / sqrt(
                                 #j10] <- tanh(atanh(toto) + qnorm(0.975) / sqrt(
                                 eno - 3))
                          CORR[jexp, jobs, 3, j3, j4, j5, j6, j7, j8, j9, 
                               j10] <- tanh(atanh(toto) + qnorm(conf_low) / sqrt(
                                 #j10] <- tanh(atanh(toto) + qnorm(0.025) / sqrt(
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
  }
  #
  dim(CORR) <- c(nexp, nobs, nvals, dimsvar[-c(posloop, poscor)])
  #  
  #  Output
  # ~~~~~~~~
  #
  CORR
}

.Corr <- function(exp, obs, siglev = 0.95, method = 'pearson', 
                  conf = TRUE, pval = TRUE) {

  # Check 'method'
  if (!(method %in% c("kendall", "spearman", "pearson"))) {
    stop("Parameter 'method' must be one of 'kendall', 'spearman' or 'pearson'.")
    #  Check 'siglev'
    if (siglev > 1 || siglev < 0) {
      stop("Parameter 'siglev' must be higher than O and lower than 1.")
    }
  }

  p_val <- NULL
  conflow <- NULL
  confhigh <- NULL
  ens_mean <- rowMeans(exp)
  CORR <- cor(obs, ens_mean, use = "pairwise.complete.obs", method = method)
  if (pval || conf) {
    if (method == "kendall" || method == "spearman") {
      eno <- Eno(rank(obs), 1)                          
    } else if (method == "pearson") {
      eno <- Eno(obs, 1)                          
    }
  }
  if (pval && (method == "pearson")) {
    t <- CORR * sqrt((eno - 2) / (1 - (CORR ^ 2)))
    p_val <- 1 - pt(t, eno - 2)
  } 
  if (conf && (method == "pearson")) {
    conf_low <- (1 - siglev) / 2
    conf_high <- 1 - conf_low
    conf_int <- c(tanh(atanh(CORR) + qnorm(conf_low) / sqrt(eno - 3)), 
                  tanh(atanh(CORR) + qnorm(conf_high) / sqrt(eno - 3)))
    conf_int <- conf_int[!is.na(CORR)]
    conflow <- conf_int[1]
    confhigh <- conf_int[2]
  } 
  #  Output
  invisible(result <- list(corr = CORR, p_val = p_val, conf_low = conflow, conf_high = confhigh))
}
