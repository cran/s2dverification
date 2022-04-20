#'Computes Root Mean Square Error
#'
#'Computes the root mean square error for an array of forecasts, var_exp and 
#'an array of observations, var_obs, which should have the same dimensions 
#'except along the posloop dimension where the lengths can be different, with 
#'the number of experiments/models for var_exp (nexp) and the number of 
#'obserational datasets for var_obs (nobs).\cr
#'The RMSE is computed along the posRMS dimension which should correspond to 
#'the startdate dimension.\cr
#'If compROW is given, the RMSE is computed only if rows along the compROW 
#'dimension are complete between limits[1] and limits[2], i.e. there are no 
#'NAs between limits[1] and limits[2]. This option can be activated if the 
#'user wishes to account only for the forecasts for which observations are 
#'available at all leadtimes.\cr
#'Default: limits[1] = 1 and limits[2] = length(compROW dimension).\cr
#'The confidence interval relies on a chi2 distribution.\cr\cr
#'.RMS provides the same functionality but taking a matrix of ensemble 
#'members as input (exp).
#'
#'@param var_exp Matrix of experimental data.
#'@param var_obs Matrix of observational data, same dimensions as var_exp 
#'  except along posloop dimension, where the length can be nobs instead of 
#'  nexp.
#'@param posloop Dimension nobs and nexp.
#'@param posRMS Dimension along which RMSE are to be computed (the dimension 
#'  of the start dates).
#'@param compROW Data taken into account only if (compROW)th row is complete.\cr
#'  Default = NULL.
#'@param limits Complete between limits[1] & limits[2]. Default = NULL.
#'@param siglev Confidence level of the computed confidence interval. 0.95 
#'  by default.
#'@param conf Whether to compute confidence interval or not. TRUE by default.
#'@param exp N by M matrix of N forecasts from M ensemble members.
#'@param obs Vector of the corresponding observations of length N.
#'
#'@return 
#'RMS: Array with dimensions:\cr
#'c(length(posloop) in var_exp, length(posloop) in var_obs, 1 or 3, all 
#'  other dimensions of var_exp & var_obs except posRMS).\cr
#'The 3rd dimension corresponds to the lower limit of the  95\% confidence 
#'  interval (only present if \code{conf = TRUE}), the RMSE, and the upper 
#'  limit of the 95\% confidence interval (only present if 
#'  \code{conf = TRUE}).\cr\cr
#'.RMS: 
#'\item{$rms}{
#'The root mean square error, 
#'}
#'\item{$conf_low}{
#'  Corresponding to the lower limit of the \code{siglev}\% confidence interval 
#'  (only present if \code{conf = TRUE}) for the rms.
#'} 
#'\item{$conf_high}{
#'  Corresponding to the upper limit of the \code{siglev}\% confidence interval 
#'  (only present if \code{conf = TRUE}) for the rms.
#'} 
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-05  (V. Guemas)  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens)  -  Formatting to R CRAN\cr
#'1.1  -  2017-02  (A. Hunter)  -  Adapted to veriApply()
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
#'clim <- Clim(sampleData$mod, sampleData$obs)
#'ano_exp <- Ano(sampleData$mod, clim$clim_exp)
#'ano_obs <- Ano(sampleData$obs, clim$clim_obs)
#'runmean_months <- 12
#'dim_to_smooth <- 4  # Smooth along lead-times
#'smooth_ano_exp <- Smoothing(ano_exp, runmean_months, dim_to_smooth)
#'smooth_ano_obs <- Smoothing(ano_obs, runmean_months, dim_to_smooth)
#'dim_to_mean <- 2  # Mean along members
#'# Discard start-dates for which some leadtimes are missing
#'required_complete_row <- 3
#'leadtimes_per_startdate <- 60
#'rms <- RMS(Mean1Dim(smooth_ano_exp, dim_to_mean), 
#'           Mean1Dim(smooth_ano_obs, dim_to_mean), 
#'           compROW = required_complete_row, 
#'           limits = c(ceiling((runmean_months + 1) / 2), 
#'                      leadtimes_per_startdate - floor(runmean_months / 2)))
#'  \donttest{
#'PlotVsLTime(rms, toptitle = "Root Mean Square Error", ytitle = "K", 
#'            monini = 11, limits = NULL, listexp = c('CMIP5 IC3'), 
#'            listobs = c('ERSST'), biglab = FALSE, hlines = c(0), 
#'            fileout = 'tos_rms.eps')
#'  }
#'# The following example uses veriApply combined with .RMS instead of RMS
#'  \dontrun{
#'require(easyVerification)
#'RMS2 <- s2dverification:::.RMS
#'rms2 <- veriApply("RMS2",
#'                  smooth_ano_exp,
#'                  # see ?veriApply for how to use the 'parallel' option
#'                  Mean1Dim(smooth_ano_obs, dim_to_mean),
#'                  tdim = 3, ensdim = 2)
#'  }
#'
#'@rdname RMS
#'@export
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
#'@rdname RMS
#'@importFrom stats qchisq
#'@export
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
