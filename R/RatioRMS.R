#'Computes the Ratio Between The RMSE of Two Experiments
#'
#'Calculates the ratio of the RMSE for two forecasts of the same observations.\cr
#'The ratio RMSE(ens, obs) / RMSE(ens.ref, obs) is output.\cr
#'The p-value is provided by a two-sided Fischer test.\cr\cr
#'.RatioRMS provides the same functionality but taking two matrices of 
#'ensemble members (ens and ens.ref) as input.
#'
#'@param var_exp1 Array of experimental data 1.
#'@param var_exp2 Array of experimental data 2.
#'@param var_obs Array of observations.
#'@param posRMS Dimension along which the RMSE are to be computed = the 
#'  position of the start dates.
#'@param pval Whether to compute the p-value of Ho : RMSE1/RMSE2 = 1 or not. 
#'  TRUE by default.
#'@param exp Matrix of experimental data 1.
#'@param exp_ref Matrix of experimental data 2.
#'@param obs Vector of observations.
#'
#'@return RatioRMS:\cr
#'Matrix with the same dimensions as var_exp1/var_exp2/var_obs except along 
#'  posRMS where the dimension has length 2 if 'pval = TRUE', or 1 otherwise. 
#'  The dimension of length 2 corresponds to the ratio between the RMSE 
#'  (RMSE1/RMSE2) and the p-value of the two-sided Fisher test with Ho: 
#'  RMSE1/RMSE2 = 1.\cr\cr
#'.RatioRMS:\cr
#'  \itemize{
#'    \item{ratiorms}{The ratio of the RMSE of the two experimental datasets}
#'    \item{p_val}{The p-value}
#'  }
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-11  (V. Guemas)  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens)  -  Formatting to R CRAN\cr
#'1.1  -  2017-02  (A. Hunter)  -  Adapted to veriApply()
#'@examples
#'# See examples on Load() to understand the first lines in this example
#'  \dontrun{
#'data_path <- system.file('sample_data', package = 's2dverification')
#'expA <- list(name = 'experiment', path = file.path(data_path,
#'             'model/$EXP_NAME$/$STORE_FREQ$_mean/$VAR_NAME$_3hourly',
#'             '$VAR_NAME$_$START_DATE$.nc'))
#'obsX <- list(name = 'observation', path = file.path(data_path,
#'             '$OBS_NAME$/$STORE_FREQ$_mean/$VAR_NAME$',
#'             '$VAR_NAME$_$YEAR$$MONTH$.nc'))
#'
#'# Now we are ready to use Load().
#'startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
#'sampleData <- Load('tos', list(expA), list(obsX), startDates,
#'                   output = 'lonlat', latmin = 27, latmax = 48, 
#'                   lonmin = -12, lonmax = 40)
#'  }
#'  \dontshow{
#'startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
#'sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'),
#'                                                c('observation'), startDates,
#'                                                output = 'lonlat',
#'                                                latmin = 27, latmax = 48,
#'                                                lonmin = -12, lonmax = 40)
#'  }
#'# Compute DJF seasonal means and anomalies.
#'leadtimes_dimension <- 4
#'initial_month <- 11
#'mean_start_month <- 12
#'mean_stop_month <- 2                                
#'sampleData$mod <- Season(sampleData$mod, leadtimes_dimension, initial_month, 
#'                         mean_start_month, mean_stop_month)
#'sampleData$obs <- Season(sampleData$obs, leadtimes_dimension, initial_month, 
#'                         mean_start_month, mean_stop_month)                              
#'clim <- Clim(sampleData$mod, sampleData$obs)
#'ano_exp <- Ano(sampleData$mod, clim$clim_exp)
#'ano_obs <- Ano(sampleData$obs, clim$clim_obs)
#'# Generate two experiments with 2 and 1 members from the only experiment 
#'# available in the sample data. Take only data values for a single forecast
#'# time step.
#'ano_exp_1 <- Subset(ano_exp, 'member', c(1, 2))
#'ano_exp_2 <- Subset(ano_exp, 'member', c(3))
#'ano_exp_1 <- Subset(ano_exp_1, c('dataset', 'ftime'), list(1, 1), drop = 'selected')
#'ano_exp_2 <- Subset(ano_exp_2, c('dataset', 'ftime'), list(1, 1), drop = 'selected')
#'ano_obs <- Subset(ano_obs, c('dataset', 'ftime'), list(1, 1), drop = 'selected')
#'# Compute ensemble mean and provide as inputs to RatioRMS.
#'rrms <- RatioRMS(Mean1Dim(ano_exp_1, 1), 
#'                 Mean1Dim(ano_exp_2, 1), 
#'                 Mean1Dim(ano_obs, 1))
#'# Plot the RatioRMS for the first forecast time step.
#'PlotEquiMap(rrms[1, , ], sampleData$lon, sampleData$lat, 
#'            toptitle = 'Ratio RMSE')
#'
#'# The following example uses veriApply combined with .RatioRMS instead of RatioRMS
#'  \dontrun{
#'require(easyVerification)
#'# The name of the function has to end in 'ss' in order for veriApply() to 
#'# detect it as a skill score.
#'RatioRMSss <- s2dverification:::.RatioRMS
#'rrms2 <- veriApply("RatioRMSss", ano_exp_1,
#'                   # see ?veriApply for how to use the 'parallel' option
#'                   Mean1Dim(ano_obs, 1),
#'                   ano_exp_2,
#'                   tdim = 2, ensdim = 1)
#'  }
#'@rdname RatioRMS
#'@export
RatioRMS <- function(var_exp1, var_exp2, var_obs, posRMS = 1, pval = TRUE) {
  #
  #  Enlarge var_exps & var_obs to 10 dim + move posRMS to 1st pos 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimsvar <- dim(var_exp1)
  for (iind in 1:length(dimsvar)) {
    if (dim(var_exp2)[iind] != dimsvar[iind] | 
        dim(var_obs)[iind] != dimsvar[iind]) {
      stop("all input vars should have the same dimensions")
    }
  } 
  enlvarexp1 <- Enlarge(var_exp1, 10)
  enlvarexp2 <- Enlarge(var_exp2, 10)
  enlvarobs <- Enlarge(var_obs, 10)
  posaperm <- 1:10
  posaperm[1] <- posRMS
  posaperm[posRMS] <- 1
  permvarexp1 <- aperm(enlvarexp1, posaperm)
  permvarexp2 <- aperm(enlvarexp2, posaperm)
  permvarobs <- aperm(enlvarobs, posaperm)
  dimsaperm <- dim(permvarexp1)
  #
  #  RMS ratio and its pvalue computation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (pval) {
    nvals <- 2
  } else {
    nvals <- 1
  }
  enlratiorms <- array(dim = c(nvals, dimsaperm[2:10]))
  dif1 <- permvarexp1 - permvarobs
  dif2 <- permvarexp2 - permvarobs
  rms1 <- Mean1Dim(dif1 ** 2, 1, narm = TRUE) ** 0.5
  rms2 <- Mean1Dim(dif2 ** 2, 1, narm = TRUE) ** 0.5
  rms2[which(abs(rms2) <= (max(abs(rms2), na.rm = TRUE) / 1000))] <- max(
    abs(rms2), na.rm = TRUE) / 1000
  enlratiorms[1, , , , , , , , , ] <- (rms1 / rms2)
  if (pval) {
    eno1 <- Eno(dif1, 1)
    eno2 <- Eno(dif2, 1)
    F <- (eno1 * (rms1) ** 2 / (eno1 - 1)) / (eno2 * (rms2) ** 2 / (eno2 - 1))
    F[which(F < 1)] <- 1 / F[which(F < 1)]
    for (j2 in 1:dimsaperm[2]) {
      for (j3 in 1:dimsaperm[3]) {
        for (j4 in 1:dimsaperm[4]) {
          for (j5 in 1:dimsaperm[5]) {
            for (j6 in 1:dimsaperm[6]) {
              for (j7 in 1:dimsaperm[7]) {
                for (j8 in 1:dimsaperm[8]) {
                  for (j9 in 1:dimsaperm[9]) {
                    for (j10 in 1:dimsaperm[10]) {
                      l1 <- eno1[j2, j3, j4, j5, j6, j7, j8, j9, j10]
                      l2 <- eno2[j2, j3, j4, j5, j6, j7, j8, j9, j10]
                      if (!is.na(l1) && !is.na(l2) && l1 > 2 && l2 > 2) {
                        enlratiorms[2, j2, j3, j4, j5, j6, j7, j8, j9, 
                                    j10] <- (1 - pf(F[j2, j3, j4, j5, j6, j7, j8, j9, j10], 
                                                    l1 - 1, l2 - 1)) * 2
                      } else {
                        enlratiorms[1, j2, j3, j4, j5, j6, j7, j8, j9, j10] <- NA
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
  
  enlratiorms <- aperm(enlratiorms, posaperm)
  if (pval) {
    dimsvar[posRMS] <- 2
  } else {
    dimsvar[posRMS] <- 1
  }
  dim(enlratiorms) <- dimsvar
  #  
  #  Output
  # ~~~~~~~~
  #
  enlratiorms
}
#'@rdname RatioRMS
#'@export
.RatioRMS <- function(exp, exp_ref, obs, pval = TRUE) {
  #
  #  RMS ratio and its pvalue computation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dif1 <- rowMeans(exp, na.rm = TRUE) - obs
  dif2 <- rowMeans(exp_ref, na.rm = TRUE) - obs
  rms1 <- mean(dif1 ** 2, na.rm = TRUE) ** 0.5
  rms2 <- mean(dif2 ** 2, na.rm = TRUE) ** 0.5
  rms2[which(abs(rms2) <= (max(abs(rms2), na.rm = TRUE) / 1000))] <- max(
                           abs(rms2), na.rm = TRUE) / 1000
  enlratiorms <- (rms1 / rms2)
  if (pval) {
    eno1 <- Eno(dif1, 1)
    eno2 <- Eno(dif2, 1)
    F <- (eno1 * (rms1) ** 2 / (eno1 - 1)) / (eno2 * (rms2) ** 2 / (eno2 - 1))
    F[which(F < 1)] <- 1 / F[which(F < 1)]
    if (!is.na(eno1) && !is.na(eno2) && eno1 > 2 && eno2 > 2) {
      p_val <- (1 - pf(F, eno1 - 1, eno2 - 1)) * 2
    } 
  }
  
  # Output
  list(ratiorms = enlratiorms, p_val = p_val)
}
