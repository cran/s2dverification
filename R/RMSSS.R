#'Computes Root Mean Square Skill Score
#'
#'Computes the root mean square error skill score between an array of 
#'forecasts, var_exp and an array of observations, var_obs, which should 
#'have the same dimensions except along posloop where the lengths can be 
#'different, with the number of experiments/models for var_exp (nexp) and 
#'the number of obserational datasets for var_obs (nobs).\cr
#'RMSSS computes the Root Mean Square Skill Score of each jexp in 1:nexp 
#'against each jobs in 1:nobs which gives nexp x nobs RMSSS for each other 
#'grid point of the matrix (each latitude/longitude/level/leadtime).\cr
#'The RMSSS are computed along the posRMS dimension which should correspond 
#'to the startdate dimension.\cr
#'The p-value is optionally provided by a one-sided Fisher test.\cr\cr
#'.RMSSS provides the same functionality but taking a matrix of ensemble 
#'members as input (exp).
#'
#'@param var_exp Array of experimental data.
#'@param var_obs Array of observational data, same dimensions as var_exp 
#'  except along posloop dimension, where the length can be nobs instead of 
#'  nexp.
#'@param posloop Dimension nobs and nexp.
#'@param posRMS Dimension along which the RMSE are to be computed (the 
#'  dimension of the start dates).
#'@param pval Whether to compute or not the p-value of the test Ho : 
#'  RMSSS = 0. TRUE by default.
#'@param exp N by M matrix of N forecasts from M ensemble members.
#'@param obs Vector of the corresponding observations of length N. 
#'
#'@return RMSSS: Array with dimensions:\cr
#'  c(length(posloop) in var_exp, length(posloop) in var_obs, 1 or 2, 
#'  all other dimensions of var_exp & var_obs except posRMS).\cr
#'The 3rd dimension corresponds to the RMSSS and, if \code{pval = TRUE}, 
#'  the p-value of the one-sided Fisher test with Ho: RMSSS = 0.\cr\cr
#'.RMSSS:
#'  \itemize{
#'    \item{$rmsss}{
#'    The RMSSS.
#'    }
#'    \item{$p_val}{
#'    Corresponds to the p values (only present if \code{pval = TRUE}) for 
#'    the RMSSS.
#'    }
#'  }
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2012-04  (V. Guemas)  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens)  -  Formatting to R CRAN\cr
#'1.1  -  2017-02  (A. Hunter)  -  Adapted to veriApply()
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
#'clim <- Clim(sampleData$mod, sampleData$obs)
#'ano_exp <- Ano(sampleData$mod, clim$clim_exp)
#'ano_obs <- Ano(sampleData$obs, clim$clim_obs)
#'rmsss <- RMSSS(Mean1Dim(ano_exp, 2), Mean1Dim(ano_obs, 2))
#'rmsss_plot <- array(dim = c(dim(rmsss)[1:2], 4, dim(rmsss)[4]))
#'rmsss_plot[, , 2, ] <- rmsss[, , 1, ]
#'rmsss_plot[, , 4, ] <- rmsss[, , 2, ]
#'  \donttest{
#'PlotVsLTime(rmsss_plot, toptitle = "Root Mean Square Skill Score", ytitle = "", 
#'            monini = 11, limits = c(-1, 1.3), listexp = c('CMIP5 IC3'), 
#'            listobs = c('ERSST'), biglab = FALSE, hlines = c(-1, 0, 1), 
#'            fileout = 'tos_rmsss.eps')
#'  }
#'# The following example uses veriApply combined with .RMSSS instead of RMSSS
#'  \dontrun{
#'require(easyVerification)
#'RMSSS2 <- s2dverification:::.RMSSS
#'rmsss2 <- veriApply("RMSSS2", ano_exp,
#'                    # see ?veriApply for how to use the 'parallel' option
#'                    Mean1Dim(ano_obs, 2),
#'                    tdim = 3, ensdim = 2)
#'  }
#'@rdname RMSSS
#'@export
RMSSS <- function(var_exp, var_obs, posloop = 1, posRMS = 2, pval = TRUE) {
  #
  #  Enlarge var_exp & var_obs & clim to 10 dim + move posloop & posRMS  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimsvar <- dim(var_exp)
  for (iind in 1:length(dimsvar)) {
    if (iind !=posloop & dim(var_obs)[iind] != dimsvar[iind]) { 
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
  #  RMSSS and its pvalue computation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (pval) {
    nvals <- 2
  } else {
    nvals <- 1
  }
  enlRMSSS <- array(dim = c(nexp, nobs, nvals, dimsaperm[3:10]))
  
  for (jexp in 1:nexp) {
    for (jobs in 1:nobs) {
      dif1 <- array(dim = dimsaperm[-1])
      dif2 <- array(dim = dimsaperm[-1])
      dif1[, , , , , , , , ] <- permvarexp[jexp, , , , , , , , 
                                           , ] - permvarobs[jobs, , , , , , , , , ]
      dif2[, , , , , , , , ] <- permvarobs[jobs, , , , , , , , , ]
      rms1 <- Mean1Dim(dif1 ** 2, 1, narm = TRUE) ** 0.5
      rms2 <- Mean1Dim(dif2 ** 2, 1, narm = TRUE) ** 0.5
      rms2[which(abs(rms2) <= (max(abs(rms2), na.rm = TRUE) / 1000))] <- max(abs(
        rms2), na.rm = TRUE) / 1000
      enlRMSSS[jexp, jobs, 1, , , , , , , , ] <- 1 - (rms1 / rms2)
      if (pval) {
        eno1 <- Eno(dif1, 1)
        eno2 <- Eno(dif2, 1)
        F <- (eno2 * (rms2) ** 2 / (eno2 - 1)) / (eno1 * (rms1) ** 2 / (eno1 - 1))
        for (j3 in 1:dimsaperm[3]) {
          for (j4 in 1:dimsaperm[4]) {
            for (j5 in 1:dimsaperm[5]) {
              for (j6 in 1:dimsaperm[6]) {
                for (j7 in 1:dimsaperm[7]) {
                  for (j8 in 1:dimsaperm[8]) {
                    for (j9 in 1:dimsaperm[9]) {
                      for (j10 in 1:dimsaperm[10]) {
                        l1 <- eno1[j3, j4, j5, j6, j7, j8, j9, j10]
                        l2 <- eno2[j3, j4, j5, j6, j7, j8, j9, j10]
                        if (is.na(l1) == FALSE & is.na(l2) == FALSE & l1 > 2 & l2 > 2) {
                          enlRMSSS[jexp, jobs, 2, j3, j4, j5, j6, j7, j8, j9, 
                                   j10] <- 1 - pf(F[j3, j4, j5, j6, j7, j8, j9, 
                                                    j10], l1 - 1, l2 - 1)
                        } else {
                          enlRMSSS[jexp, jobs, 1, j3, j4, j5, j6, j7, j8, j9,
                                   j10] <- NA
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
  
  dim(enlRMSSS) <- c(nexp, nobs, nvals, dimsvar[-c(posloop, posRMS)])
  #
  #  Output
  # ~~~~~~~~
  #
  enlRMSSS
}
#'@rdname RMSSS
#'@importFrom stats pf
#'@export
.RMSSS <- function(exp, obs, pval = TRUE) {
  dif2 <- obs
  dif1 <- rowMeans(exp) - obs
  rms1 <- mean(dif1 ** 2, na.rm = TRUE) ** 0.5
  rms2 <- mean(dif2 ** 2, na.rm = TRUE) ** 0.5
  rms2[abs(rms2) <= abs(rms2) / 1000] <- abs(rms2) / 1000
  rmsss <- c(1 - (rms1 / rms2))
  if (pval == TRUE) {
    eno1 <- Eno(dif1, 1)
    eno2 <- Eno(dif2, 1)
    F.stat <- (eno2 * (rms2) ** 2 / (eno2 - 1)) / (eno1 * (rms1) ** 2 / (eno1 - 1))
    if (is.na(eno1) == FALSE && is.na(eno2) == FALSE && eno1 > 2 && eno2 > 2) {
      p_val <- 1 - pf(F.stat, eno1 - 1, eno2 - 1)
    }
  } else {
    p_val <- NA
  }

  #
  #  Output
  # ~~~~~~~~
  #
  list(rmsss = rmsss, p_val = p_val)
}
