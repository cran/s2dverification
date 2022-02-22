#'Computes InterQuartile Range, Maximum-Minimum, Standard Deviation and 
#'Median Absolute Deviation of the Ensemble Members
#'
#'Computes the InterQuartile Range, the Maximum minus Mininum, the Standard 
#'Deviation and the Median Absolute Deviation along the list of dimensions 
#'provided by the posdim argument (typically along the ensemble member and 
#'start date dimension).\cr
#'The confidence interval is optionally computed by bootstrapping.
#'
#'@param var Matrix of any number of dimensions up to 10.
#'@param posdim List of dimensions along which to compute IQR/MaxMin/SD/MAD.
#'@param narm TRUE/FALSE if NA removed/kept for computation. Default = TRUE.
#'@param siglev Confidence level of the computed confidence interval. 
#'  0.95 by default.
#'@param conf Whether to compute the confidence intervals or not. 
#'  TRUE by default.
#'
#'@details
#'Example:\cr
#'--------\cr
#'To compute IQR, Max-Min, SD & MAD accross the members and start dates of 
#'var output from \code{Load()} or \code{Ano()} or \code{Ano_CrossValid()}, 
#'call:\cr
#'  spread(var, posdim = c(2, 3), narm = TRUE)
#'
#'@return 
#'Matrix with the same dimensions as var except along the first posdim 
#'dimension which is replaced by a length 1 or 3 dimension, corresponding to 
#'the lower limit of the \code{siglev}\% confidence interval 
#'(only present if \code{conf = TRUE}), the spread, and the upper limit of 
#'the \code{siglev}\% confidence interval (only present if \code{conf = TRUE}) 
#'for each experiment/leadtime/latitude/longitude.
#'  \item{$iqr}{
#'    InterQuartile Range.
#'  }
#'  \item{$maxmin}{
#'    Maximum - Minimum.
#'  }
#'  \item{$sd}{
#'    Standard Deviation.
#'  }
#'  \item{$mad}{
#'    Median Absolute Deviation.
#'  } 
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-03  (V. Guemas, \email{virginie.guemas@@bsc.es})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@@bsc.es})  -  Formatting to CRAN
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
#'clim <- Clim(sampleData$mod, sampleData$obs)
#'ano_exp <- Ano(sampleData$mod, clim$clim_exp)
#'runmean_months <- 12
#'dim_to_smooth <- 4  # Smooth along lead-times
#'smooth_ano_exp <- Smoothing(ano_exp, runmean_months, dim_to_smooth)
#'smooth_ano_exp_m_sub <- smooth_ano_exp - InsertDim(Mean1Dim(smooth_ano_exp, 2, 
#'                        narm = TRUE), 2, dim(smooth_ano_exp)[2])
#'spread <- Spread(smooth_ano_exp_m_sub, c(2, 3))
#'  \donttest{
#'PlotVsLTime(spread$iqr, 
#'            toptitle = "Inter-Quartile Range between ensemble members",
#'            ytitle = "K", monini = 11, limits = NULL, 
#'            listexp = c('CMIP5 IC3'), listobs = c('ERSST'), biglab = FALSE, 
#'            hlines = c(0), fileout = 'tos_iqr.eps')
#'PlotVsLTime(spread$maxmin, toptitle = "Maximum minus minimum of the members", 
#'            ytitle = "K", monini = 11, limits = NULL, 
#'            listexp = c('CMIP5 IC3'), listobs = c('ERSST'), biglab = FALSE, 
#'            hlines = c(0), fileout = 'tos_maxmin.eps')
#'PlotVsLTime(spread$sd, toptitle = "Standard deviation of the members", 
#'            ytitle = "K", monini = 11, limits = NULL, 
#'            listexp = c('CMIP5 IC3'), listobs = c('ERSST'), biglab = FALSE, 
#'            hlines = c(0), fileout = 'tos_sd.eps')
#'PlotVsLTime(spread$mad, toptitle = "Median Absolute Deviation of the members",
#'            ytitle = "K", monini = 11, limits = NULL, 
#'            listexp = c('CMIP5 IC3'), listobs = c('ERSST'), biglab = FALSE, 
#'            hlines = c(0), fileout = 'tos_mad.eps')
#'  }
#'
#'@importFrom stats IQR sd mad runif quantile
#'@export
Spread <- function(var, posdim = 2, narm = TRUE, siglev = 0.95, conf = TRUE) {
  # 
  #  Enlarge the size of var to 10 and move all posdim to first position 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  dimsvar <- dim(var)
  enlvar <- Enlarge(var, 10)
  outdim <- c(dimsvar, array(1, dim = (10 - length(dimsvar))))
  posaperm <- 1:10
  posaperm[(length(posdim) + 1):10] <- posaperm[-posdim]
  posaperm[1:length(posdim)] <- posdim
  enlvar <- aperm(enlvar, posaperm)
  dimsaperm <- outdim[posaperm]
  outdim <- 1
  for (jpos in 1:length(posdim)) {
    outdim <- outdim * dimsaperm[jpos]
  }
  outdim <- c(outdim, dimsaperm[(length(posdim) + 1):10])
  enlvar <- array(enlvar, dim = outdim)
  enlvar <- Enlarge(enlvar, 10)
  outdim <- c(outdim, array(1, dim = (10 - length(outdim))))
  if (conf) {
    nvals <- 3
    posvar <- 2
  } else {
    nvals <- 1
    posvar <- 1
  }
  outdim[1] <- nvals
  # 
  #  IQR / Max-Min / STD / Median computation
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  iqrvar <- array(dim = outdim)
  maxminvar <- array(dim = outdim)
  sdvar <- array(dim = outdim)
  medvar <- array(dim = outdim)
  if (conf) {
    conf_low <- (1 - siglev) / 2
    conf_high <- 1 - conf_low
  }
  for (j2 in 1:outdim[2]) {
    for (j3 in 1:outdim[3]) {
      for (j4 in 1:outdim[4]) {
        for (j5 in 1:outdim[5]) {
          for (j6 in 1:outdim[6]) {
            for (j7 in 1:outdim[7]) {
              for (j8 in 1:outdim[8]) { 
                for (j9 in 1:outdim[9]) {
                  for (j10 in 1:outdim[10]) {
                    # 
                    #  Compute spread along the posdim dimension 
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    # 
                    whichspread <- enlvar[, j2, j3, j4, j5, j6, j7, j8, j9, j10]
                    iqrvar[posvar, j2, j3, j4, j5, j6, j7, j8, j9, 
                           j10] <- IQR(whichspread, na.rm = narm)
                    maxminvar[posvar, j2, j3, j4, j5, j6, j7, j8, j9, 
                              j10] <- max(whichspread, 
                                      na.rm = narm) - min(whichspread, 
                                                          na.rm = narm)
                    sdvar[posvar, j2, j3, j4, j5, j6, j7, j8, j9, 
                          j10] <- sd(whichspread, na.rm = narm)
                    medvar[posvar, j2, j3, j4, j5, j6, j7, j8, j9, 
                           j10] <- mad(whichspread, na.rm = narm)           
                    if (conf) {
                      # 
                      #  Bootstrapping and recomputing the spread
                      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      # 
                      iqr_bs <- c()
                      maxmin_bs <- c()
                      sd_bs <- c()
                      mad_bs <- c()
                      Num <- length(whichspread)
                      for (jmix in 1:100) {
                        drawings <- round(runif(Num, 0.5, Num + 0.5))
                        iqr_bs <- c(iqr_bs, IQR(whichspread[drawings],  
                                    na.rm = narm))
                        maxmin_bs <- c(maxmin_bs, max(whichspread[drawings],
                                       na.rm=narm) - min(whichspread[drawings],
                                                         na.rm = narm))
                        sd_bs <- c(sd_bs, sd(whichspread[drawings], na.rm = narm))
                        mad_bs <- c(mad_bs, mad(whichspread[drawings], 
                                                na.rm = narm))
                      }
                      #
                      #  95% confidence interval obtained from the bootstrapping 
                      #  results
                      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      #  (This confidence interval is itself quite uncertain but 
                      #   it is not possible to produce miracle with 5 values ...)
                      #
                      iqrvar[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- quantile(iqr_bs, conf_low, na.rm = narm)
                      iqrvar[3, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- quantile(iqr_bs, conf_high, na.rm = narm)
                      maxminvar[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- maxminvar[2, j2, j3, j4, j5, j6, j7, j8, 
                                        j9, j10] + (quantile(maxmin_bs, conf_low, 
                                        na.rm = narm) - quantile(maxmin_bs, 
                                        conf_high, na.rm = narm)) / 2
                      maxminvar[3, j2, j3, j4, j5, j6, j7, j8, j9, 
                                j10] <- maxminvar[2, j2, j3, j4, j5, j6, j7, j8, 
                                        j9, j10] - (quantile(maxmin_bs, conf_low, 
                                        na.rm = narm) - quantile(maxmin_bs, 
                                        conf_high, na.rm = narm)) / 2
                      sdvar[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                            j10] <- quantile(sd_bs, conf_low, na.rm = narm)
                      sdvar[3, j2, j3, j4, j5, j6, j7, j8, j9, 
                            j10] <- quantile(sd_bs, conf_high, na.rm = narm)
                      medvar[1, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- medvar[2, j2, j3, j4, j5, j6, j7, j8, j9, 
                                     j10] + (quantile(mad_bs, conf_low, 
                                     na.rm = narm) - quantile(mad_bs, conf_high, 
                                     na.rm = narm))
                      medvar[3, j2, j3, j4, j5, j6, j7, j8, j9, 
                             j10] <- medvar[2, j2, j3, j4, j5, j6, j7, j8, j9, 
                                     j10] - (quantile(mad_bs, conf_low, 
                                     na.rm = narm) - quantile(mad_bs, conf_high, 
                                     na.rm = narm))
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
  posaperm <- posdim[1]
  dimsvar <- dimsvar[-posdim]
  if (posdim[1] > 1) {
    posaperm <- c(posaperm, 1:(posdim[1] - 1)) 
  }
  if (posdim[1] <= length(dimsvar)) {
    posaperm <- c(posaperm, (posdim[1] + 1):(length(posdim) + length(dimsvar)))
  }
  while (length(posaperm) < 10) {
    posaperm <- c(posaperm,length(posaperm) + 1) 
  }
  outdim <- outdim[posaperm]
  outdim <- outdim[1:(length(dimsvar) + 1)]
  iqrvar <- aperm(iqrvar, posaperm)
  dim(iqrvar) <- outdim
  maxminvar <- aperm(maxminvar, posaperm)
  dim(maxminvar) <- outdim
  maxminvar[which(is.infinite(maxminvar))] <- NA
  sdvar <- aperm(sdvar, posaperm)
  dim(sdvar) <- outdim
  medvar <- aperm(medvar, posaperm)
  dim(medvar) <- outdim
  medvar[which(is.infinite(medvar))] <- NA
  # 
  #  Outputs
  # ~~~~~~~~~
  # 
  invisible(list(iqr = iqrvar, maxmin = maxminvar, sd = sdvar, mad = medvar))
}
