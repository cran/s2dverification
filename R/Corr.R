#'Computes the correlation coefficient between an array of forecasts and their corresponding observations
#'
#'Calculates the correlation coefficient (Pearson, Kendall or Spearman) for 
#'an array of forecasts and observations. The input should be an array with 
#'dimensions c(no. of datasets, no. of start dates, no. of forecast times, 
#'no. of lons, no. of lats.), where the longitude and latitude dimensions are 
#'optional. The correlations are computed along the poscor dimension which 
#'should correspond to the startdate dimension. If compROW is given, the 
#'correlations are computed only if rows along the compROW dimension are 
#'complete between limits[1] and limits[2], i.e. there are no NAs between 
#'limits[1] and limits[2]. This option can be activated if the user wishes to 
#'account only for the forecasts for which observations are available at all 
#'leadtimes. \cr 
#'Default: limits[1] = 1 and limits[2] = length(compROW dimension).\cr 
#'The confidence interval is computed by a Fisher transformation.\cr 
#'The significance level relies on a one-sided student-T distribution.\cr 
#'We can modifiy the treshold of the test modifying siglev (default value=0.95).\cr\cr 
#'.Corr calculates the correlation between the ensemble mean and the 
#'observations, using an N by M matrix (exp) of forecasts and a vector of 
#'observations (obs) as input.
#'
#'@param var_exp Array of experimental data.
#'@param var_obs Array of observational data, same dimensions as var_exp 
#'  except along posloop dimension, where the length can be nobs instead of nexp.
#'@param posloop Dimension nobs and nexp.
#'@param poscor Dimension along which correlation are to be computed (the 
#'  dimension of the start dates).
#'@param compROW Data taken into account only if (compROW)th row is complete.
#'  Default = NULL.
#'@param limits Complete between limits[1] & limits[2]. Default = NULL.
#'@param siglev Significance level. Default = 0.95.
#'@param method Type of correlation: 'pearson', 'spearman' or 'kendall'. 
#'  Default='pearson'
#'@param conf Whether to compute confidence intervals (default = 'TRUE') or 
#'  not (FALSE).
#'@param pval Whether to compute statistical significance p-value (default = 'TRUE') 
#'  or not (FALSE).
#'@param exp N by M matrix of N forecasts from M ensemble members.
#'@param obs Vector of the corresponding observations of length N.
#'
#'@return 
#'Corr: Array with dimensions :\cr 
#'c(# of datasets along posloop in var_exp, # of datasets along posloop in 
#'var_obs, 4, all other dimensions of var_exp & var_obs except poscor).\cr 
#'The third dimension, of length 4 maximum, contains to the lower limit of 
#'the 95\% confidence interval, the correlation, the upper limit of the 95\% 
#'confidence interval and the 95\% significance level given by a one-sided 
#'T-test. If the p-value is disabled via \code{pval = FALSE}, this dimension 
#'will be of length 3. If the confidence intervals are disabled via 
#'\code{conf = FALSE}, this dimension will be of length 2. If both are 
#'disabled, this will be of length 2. \cr\cr
#'.Corr:
#'  \itemize{
#'    \item{$corr}{The correlation statistic.}
#'    \item{$p_val}{Corresponds to the p values for the \code{siglev}\% 
#'    (only present if \code{pval = TRUE}) for the correlation.} 
#'    \item{$conf_low}{Corresponds to the upper limit of the \code{siglev}\% 
#'    (only present if \code{conf = TRUE}) for the correlation.} 
#'    \item{$conf_high}{Corresponds to the lower limit of the \code{siglev}\% 
#'    (only present if \code{conf = TRUE}) for the correlation.}
#'  }
#'
#'@keywords datagen
#'@author History:\cr
#'  0.1  -  2011-04  (V. Guemas, \email{vguemas@@bsc.es})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@@bsc.es})  -  Formatting to R CRAN\cr
#'1.1  -  2014-10  (M. Menegoz, \email{martin.menegoz@@bsc.es})  -  Adding siglev argument\cr
#'1.2  -  2015-03  (L.P. Caron, \email{louis-philippe.caron@@bsc.es})  - Adding method argument\cr
#'1.3  -  2017-02  (A. Hunter, \email{alasdair.hunter@@bsc.es})  -  Adapted to veriApply()
#'@examples
#'# Load sample data as in Load() example: 
#'example(Load) 
#'clim <- Clim(sampleData$mod, sampleData$obs) 
#'ano_exp <- Ano(sampleData$mod, clim$clim_exp) 
#'ano_obs <- Ano(sampleData$obs, clim$clim_obs) 
#'runmean_months <- 12 
#'dim_to_smooth <- 4  
#'# Smooth along lead-times   
#'smooth_ano_exp <- Smoothing(ano_exp, runmean_months, dim_to_smooth) 
#'smooth_ano_obs <- Smoothing(ano_obs, runmean_months, dim_to_smooth) 
#'dim_to_mean <- 2  # Mean along members 
#'required_complete_row <- 3  # Discard start dates which contain any NA lead-times 
#'leadtimes_per_startdate <- 60 
#'corr <- Corr(Mean1Dim(smooth_ano_exp, dim_to_mean),              
#'             Mean1Dim(smooth_ano_obs, dim_to_mean),              
#'             compROW = required_complete_row,              
#'             limits = c(ceiling((runmean_months + 1) / 2),                         
#'             leadtimes_per_startdate - floor(runmean_months / 2))) 
#'  \donttest{
#'PlotVsLTime(corr, toptitle = "correlations", ytitle = "correlation",             
#'            monini = 11, limits = c(-1, 2), listexp = c('CMIP5 IC3'),
#'            listobs = c('ERSST'), biglab = FALSE, hlines = c(-1, 0, 1),
#'            fileout = 'tos_cor.eps')
#'  }
#'
#'# The following example uses veriApply combined with .Corr instead of Corr
#'  \dontrun{
#'require(easyVerification)  
#'Corr2 <- s2dverification:::.Corr
#'corr2 <- veriApply("Corr2", 
#'                   smooth_ano_exp, 
#'                   # see ?veriApply for how to use the 'parallel' option
#'                   Mean1Dim(smooth_ano_obs, dim_to_mean), 
#'                   tdim = 3, ensdim = 2)
#'  }
#'@rdname Corr
#'@importFrom stats cor pt qnorm 
#'@export
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

#'@rdname Corr
#'@importFrom stats cor pt qnorm
#'@export
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
    t <-sqrt(CORR * CORR * (eno - 2) / (1 - (CORR ^ 2)))
    p_val <- pt(t, eno - 2, lower.tail = FALSE)
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
