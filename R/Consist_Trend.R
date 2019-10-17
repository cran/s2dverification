#'Computes Trends Using Only Model Data For Which Observations Are Available
#'
#'Computes the trend coefficients for a time series by least square fitting, 
#'together with the associated error interval for both the observational and 
#'model data.\cr
#'Provides also the detrended observational and modeled data.\cr
#'By default, the trend is computed along the second dimension of the input 
#'array, which is expected to be the start date dimension. For arrays 
#'containing multiple model members, the user will first have to calculate 
#'the ensemble average using \code{Mean1Dim()} or elsewhise (see the example).
#'
#'@param var_exp Ensemble mean of model hindcasts with dimensions:\cr
#'  c(nmod/nexp, nsdates, nltime) up to\cr
#'  c(nmod/nexp, nsdates, nltime, nlevel, nlat, nlon)
#'@param var_obs Ensemble mean of observational data with dimensions:\cr
#'  c(nobs, nsdates, nltime) up to\cr
#'  c(nobs, nsdates, nltime, nlevel, nlat, nlon)\cr
#'  Dimensions 2 to 6 should be the same as var_exp.
#'@param interval Number of months/years between 2 start dates. Default = 1. 
#'  The trends will be provided respectively in field unit per month or per year.
#'
#'@return 
#'\item{$trend}{
#'  Trend coefficients of model and observational data with dimensions:\cr
#'  c(nmod/nexp + nobs, 3, nltime) up to\cr
#'  c(nmod/nexp + nobs, 3, nltime, nlevel, nlat, nlon)\cr
#'  The length 3 dimension corresponds to the lower limit of the 95\% 
#'  confidence interval, the slope of the trends and the upper limit of the 
#'  95\% confidence interval.
#' }
#'\item{$detrendedmod}{
#'  Same dimensions as var_exp with linearly detrended values of var_exp 
#'  along the second = start date dimension.
#' }
#'\item{$detrendedobs}{
#'  Same dimensions as var_exp with linearly detrended values of var_obs 
#'  along the second = start date dimension.
#' }
#'
#'@keywords datagen
#'@author History:\cr
#'  0.1  -  2011-11  (V. Guemas, \email{vguemas@@ic3.cat})  -  Original code\cr
#'  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@@ic3.cat})  -  Formatting to R CRAN
#'@examples
#'#'# Load sample data as in Load() example:
#'example(Load)
#'clim <- Clim(sampleData$mod, sampleData$obs)
#'ano_exp <- Ano(sampleData$mod, clim$clim_exp)
#'ano_obs <- Ano(sampleData$obs, clim$clim_obs)
#'runmean_months <- 12
#'dim_to_smooth <- 4  # Smooth along lead-times
#'smooth_ano_exp <- Smoothing(ano_exp, runmean_months, dim_to_smooth)
#'smooth_ano_obs <- Smoothing(ano_obs, runmean_months, dim_to_smooth)
#'dim_to_mean <- 2  # Mean along members
#'years_between_startdates <- 5
#'trend <- Consist_Trend(Mean1Dim(smooth_ano_exp, dim_to_mean), 
#'                       Mean1Dim(smooth_ano_obs, dim_to_mean), 
#'                       years_between_startdates)
#'
#'  \donttest{
#'PlotVsLTime(trend$trend, toptitle = "trend", ytitle = "K/(5 years)", 
#'            monini = 11, limits = c(-0.8, 0.8), listexp = c('CMIP5 IC3'), 
#'            listobs = c('ERSST'), biglab = FALSE, hlines = c(0), 
#'            fileout = 'tos_consist_trend.eps')
#'PlotAno(InsertDim(trend$detrendedmod,2,1), InsertDim(trend$detrendedobs,2,1), 
#'        startDates, "Detrended tos anomalies", ytitle = 'K', 
#'        legends = 'ERSST', biglab = FALSE, fileout = 'tos_detrended_ano.eps')
#'  }
#'
#'@export
Consist_Trend <- function(var_exp, var_obs, interval = 1) {
  #
  #  Enlarge the number of dimensions of var_exp and var_obs to 6 if necessary
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimexp <- dim(var_exp)
  dimobs <- dim(var_obs)
  if (length(dimexp) < 3 | length(dimobs) < 3) {
    stop("At least 3 dim needed : c(nexp/nobs, nsdates, nltime)") 
  }
  for (jn in 2:max(length(dimexp), length(dimobs))) {  
    if (dimexp[jn] != dimobs[jn]) {
      stop("Wrong input dimensions")
    }
  }
  var_exp <- Enlarge(var_exp, 6)
  var_obs <- Enlarge(var_obs, 6)
  #
  #  Find common points to compute trends 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  tmp <- MeanListDim(var_obs, dims = 4:6, narm = TRUE)
  tmp2 <- MeanListDim(var_exp, dims = 4:6, narm = TRUE)
  nan <- Mean1Dim(tmp, 1, narm = FALSE) + Mean1Dim(tmp2, 1, narm = FALSE)
  for (jdate in 1:dimexp[2]) {
    for (jtime in 1:dimexp[3]) {
      if (is.na(nan[jdate, jtime])) {
        var_exp[, jdate, jtime, , , ] <- NA
        var_obs[, jdate, jtime, , , ] <- NA
      }
    }
  }
  #
  #  Compute trends 
  # ~~~~~~~~~~~~~~~~
  #
  trendcat <- array(dim = c(dimobs[1] + dimexp[1], 4, dim(var_obs)[3:6]))
  detrendedmod <- array(dim = dimexp)
  detrendedobs <- array(dim = dimobs)
  tmp1 <- Trend(var_exp, 2, interval = interval)
  tmp2 <- Trend(var_obs, 2, interval = interval)
  trendcat[1:dimexp[1], , , , , ] <- tmp1$trend
  trendcat[(dimexp[1] + 1):dim(trendcat)[1], , , , , ] <- tmp2$trend
  detrendedmod[] <- tmp1$detrended
  detrendedobs[] <- tmp2$detrended
  #
  #  Back to the original dimensions
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  outtrendcat <- array(dim = c(dimobs[1] + dimexp[1], 4, 
                               dimobs[3:length(dimobs)]))
  outtrendcat[] <- trendcat
  #
  #  Outputs
  # ~~~~~~~~~
  #
  invisible(list(trend = outtrendcat, detrendedmod = detrendedmod,
                 detrendedobs = detrendedobs))
}
