#'Computes the ratio between the ensemble spread and RMSE
#'
#'Arrays var_exp & var_obs should have dimensions between\cr
#'c(nmod/nexp, nmemb/nparam, nsdates, nltime)\cr
#'and\cr
#'c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)\cr
#'The ratio between the standard deviation of the members around the ensemble 
#'mean in var_exp and the RMSE between var_exp and var_obs is output for each 
#'experiment and each observational dataset.\cr
#'The p-value is provided by a one-sided Fischer test.\cr\cr
#'.RatioSDRMS provides the same functionality but taking a matrix of ensemble 
#'members as input (exp).
#'
#'@param var_exp Model data:\cr
#'  c(nmod/nexp, nmemb/nparam, nsdates, nltime) up to\cr
#'  c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)
#'@param var_obs Observational data:\cr
#'  c(nobs, nmemb, nsdates, nltime) up to\cr
#'  c(nobs, nmemb, nsdates, nltime, nlevel, nlat, nlon)
#'@param pval Whether to compute the p-value of Ho : SD/RMSE = 1 or not.
#'@param exp N by M matrix of N forecasts from M ensemble members.
#'@param obs Vector of the corresponding observations of length N.
#'
#'@return RatioSDRMS: Array with dimensions c(nexp/nmod, nobs, 1 or 2, nltime) 
#'  up to c(nexp/nmod, nobs, 1 or 2, nltime, nlevel, nlat, nlon).\cr
#'The 3rd dimension corresponds to the ratio (SD/RMSE) and the p.value 
#' (only present if \code{pval = TRUE}) of the one-sided Fisher test with 
#'Ho: SD/RMSE = 1.\cr\cr
#'.RatioSDRMS:
#'  \itemize{
#'    \item{$ratio}{
#'    The ratio of the ensemble spread and RMSE, 
#'    }
#'    \item{$p_val}{
#'    Corresponds to the p values of the ratio (only present if 
#'    \code{pval = TRUE}).
#'    }
#'  }
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-12  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau-manubens@ic3.cat})  -  Formatting to CRAN\cr
#'1.1  -  2017-02  (A. Hunter, \email{alasdair.hunter@bsc.es})  -  Adapted to veriApply()
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
#'rsdrms <- RatioSDRMS(sampleData$mod, sampleData$obs)
#'# Reorder the data in order to plot it with PlotVsLTime
#'rsdrms_plot <- array(dim = c(dim(rsdrms)[1:2], 4, dim(rsdrms)[4]))
#'rsdrms_plot[, , 2, ] <- rsdrms[, , 1, ]
#'rsdrms_plot[, , 4, ] <- rsdrms[, , 2, ]
#'  \donttest{
#'PlotVsLTime(rsdrms_plot, toptitle = "Ratio ensemble spread / RMSE", ytitle = "", 
#'            monini = 11, limits = c(-1, 1.3), listexp = c('CMIP5 IC3'), 
#'            listobs = c('ERSST'), biglab = FALSE, siglev = TRUE, 
#'            fileout = 'tos_rsdrms.eps')
#'  }
#'
#'# The following example uses veriApply combined with .RatioSDRMS instead of RatioSDRMS
#'  \dontrun{
#'require(easyVerification)  
#'RatioSDRMS2 <- s2dverification:::.RatioSDRMS
#'rsdrms2 <- veriApply("RatioSDRMS2",
#'                     sampleData$mod,
#'                     # see ?veriApply for how to use the 'parallel' option
#'                     Mean1Dim(sampleData$obs, 2),
#'                     tdim = 3, ensdim = 2)
#'  }
#'@rdname RatioSDRMS
#'@export
RatioSDRMS <- function(var_exp, var_obs, pval = TRUE) {
  #
  #  Enlarge the number of dimensions of var_exp and var_obs to 7 if necessary
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimexp <- dim(var_exp)
  dimobs <- dim(var_obs)
  if (length(dimexp) < 4 | length(dimobs) < 4) {
    stop("At least 4 dim needed : c(nexp/nobs, nmemb, nsdates, nltime)")
  }
  for (jn in 3:max(length(dimexp), length(dimobs))) {
    if (dimexp[jn] != dimobs[jn]) {
      stop("Wrong input dimensions")
    }
  }
  var_exp <- Enlarge(var_exp, 7)
  var_obs <- Enlarge(var_obs, 7)
  
  #
  #  Ratio RMSE / SD and its significance level
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  ensmeanexp <- Mean1Dim(var_exp, 2)
  ensmeanobs <- Mean1Dim(var_obs, 2)
  dimrms <- c(dimexp[1], dimobs[1], dimexp[4:length(dimexp)])
  if (pval) {
    nvals <- 2
  } else {
    nvals <- 1
  }
  dimratiormssd <- c(dimexp[1], dimobs[1], nvals, dimexp[4:length(dimexp)])
  if (length(dimrms) < 6) {
    dimrms <- c(dimrms, array(1, dim = (6 - length(dimrms))))
  }
  if (length(dimratiormssd) < 7) {
    dimenlratiormssd <- c(dimratiormssd, 
                          array(1, dim = (7 - length(dimratiormssd))))
  } else {
    dimenlratiormssd <- dimratiormssd
  }
  dif <- var_exp - InsertDim(ensmeanexp, 2, dimexp[2])
  std <- apply(array(dif, dim = c(dimexp[1], dimexp[2] * dimexp[3],
                                  dimrms[3:6])), c(1, 3, 4, 5, 6), sd, na.rm = TRUE)
  enosd <- apply(Eno(dif, 3), c(1, 3, 4, 5, 6), sum, na.rm = TRUE)
  rms <- array(dim = dimrms)
  enlratiormssd <- array(dim = dimenlratiormssd)
  for (jexp in 1:dimexp[1]) {
    for (jobs in 1:dimobs[1]) {
      dif <- ensmeanexp[jexp, , , , , ] - ensmeanobs[jobs, , , , , ]
      rms[jexp,jobs, , , , ] <- Mean1Dim(dif ** 2, 1, narm = TRUE) ** 0.5
      enorms <- array(Eno(dif, 1), dim = dimrms[3:6])
      enlratiormssd[jexp, jobs, 1, , , , ] <- std[jexp, , , , ] / rms[jexp, 
                                                                      jobs, , , , ]
      if (pval) {
        for (jltime in 1:dimrms[3]) {
          for (jlev in 1:dimrms[4]) {
            for (jlat in 1:dimrms[5]) {
              for (jlon in 1:dimrms[6]) {
                l1 <- enosd[jexp, jltime, jlev, jlat, jlon]
                l2 <- enorms[jltime, jlev, jlat, jlon]
                F <- (enosd[jexp, jltime, jlev, jlat, jlon] * (std[jexp, jltime, 
                                                                   jlev, jlat, jlon]) ** 2 / (enosd[jexp, jltime, jlev, jlat, 
                                                                                                    jlon] - 1)) / (enorms[jltime, jlev, jlat, jlon] * (rms[jexp, 
                                                                                                                                                           jobs, jltime, jlev, jlat, jlon]) ** 2 / (enorms[jltime, 
                                                                                                                                                                                                           jlev, jlat, jlon] - 1))
                if (!is.na(F) && !is.na(l1) && !is.na(l2) && l1 > 2 && l2 > 2) {
                  enlratiormssd[jexp, jobs, 2, jltime, jlev, jlat, 
                                jlon] <- 1 - pf(F, l1 - 1, l2 - 1)
                } else {
                  enlratiormssd[jexp, jobs, 1, jltime, jlev, jlat, jlon] <- NA
                }
              }
            }
          }
        }
      }
    }
  }
  dim(enlratiormssd) <- dimratiormssd
  #  
  #  Output
  # ~~~~~~~~
  #
  enlratiormssd
}
#'@rdname RatioSDRMS
#'@importFrom stats sd pf
#'@export
.RatioSDRMS <- function(exp, obs, pval = TRUE) {
  
  #
  #  Ratio RMSE / SD and its significance level
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  ensmean <- rowMeans(exp)
  enosd <- Eno(ensmean, 1)
  dif <- exp - ensmean 
  std <- sd(dif, na.rm = TRUE)
  dif <- ensmean - obs
  rms <- mean(dif ** 2, na.rm = TRUE) ** 0.5
  enorms <- Eno(dif, 1)
  enlratiormssd <- std / rms
  p_val <- 0
  if (pval) {
    l1 <- enosd
    l2 <- enorms 
    F <- (enosd * std ** 2 / (enosd - 1)) / (enorms * (rms) ** 2 / (enorms - 1))
    if (!is.na(F) && !is.na(l1) && !is.na(l2) && l1 > 2 && l2 > 2) {
      p_val <- 1 - pf(F, l1 - 1, l2 - 1)
    }
    else {
      p_val <- NA
    }
  } 
  #  
  #  Output
  # ~~~~~~~~
  #
  list(ratio = enlratiormssd, p_val = p_val)  
}
