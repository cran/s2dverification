\name{Corr}
\alias{Corr}
\title{
  Computes Correlation Skill Measure (Temporal Correlation Along Start Dates)
}
\description{
  Matrix var_exp & var_obs should have the same dimensions except along posloop   
  dimension  where the length can be different, with the number of experiments/
  models for var_exp (nexp) and the number of obserational datasets for var_obs
  (nobs).
  Corr computes the correlation skill of each jexp in 1:nexp against each 
  jobs in 1:nobs which gives nexp x nobs correlation skill measures for each 
  other grid point of the matrix (each latitude/longitude/level/leadtime). 
  The correlations are computed along the poscor dimension which should 
  correspond to the startdate dimension. If compROW is given,
  the correlations are computed only if rows along the compROW dimension 
  are complete between limits[1] and limits[2], that mean with no NA between 
  limits[1] and limits[2]. This option can be activated if the user wishes
  to account only for the forecasts for which observations are available at
  all leadtimes. 
  Default: limits[1] = 1 and limits[2] = length(compROW dimension).
  The confidence interval is computed by a Fisher transformation.
  The significance level relies on a one-sided student-T distribution.
}
\usage{
  Corr(var_exp, var_obs, posloop = 1, poscor = 2, compROW = NULL, limits = NULL)
}
\arguments{
  \item{var_exp}{
    Matrix of experimental data.
  }
  \item{var_obs}{
    Matrix of observational data, same dimensions as var_exp except along 
    posloop dimension, where the length can be nobs instead of nexp.
  }
  \item{posloop}{
    Dimension nobs and nexp.
  }
  \item{poscor}{
    Dimension along which correlation are to be computed (the dimension of 
    the start dates).
  }
  \item{compROW}{
    Data taken into account only if (compROW)th row is complete.
    Default = NULL.
  }
  \item{limits}{
    Complete between limits[1] & limits[2]. Default = NULL.
  }
}
\value{
  Matrix with dimensions :
  c(length(posloop) in var_exp, length(posloop) in var_obs, 4, all
     other dimensions of var_exp & var_obs except poscor).
  The third dimension of length 4 corresponds to the lower limit of the  95\% 
  confidence interval, the correlation, the upper limit of the 95\% 
  confidence interval and the 95\% significance level given by a one-sided 
  T-test.
}
\examples{
  startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'areave', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  clim <- Clim(sampleData$mod, sampleData$obs)
  ano_exp <- Ano(sampleData$mod, clim$clim_exp)
  ano_obs <- Ano(sampleData$obs, clim$clim_obs)
  runmean_months <- 12
  dim_to_smooth <- 4  # Smooth along lead-times
  smooth_ano_exp <- Smoothing(ano_exp, runmean_months, dim_to_smooth)
  smooth_ano_obs <- Smoothing(ano_obs, runmean_months, dim_to_smooth)
  dim_to_mean <- 2  # Mean along members
  required_complete_row <- 3  # Discard start dates which contain any NA lead-times
  leadtimes_per_startdate <- 60
  corr <- Corr(Mean1Dim(smooth_ano_exp, dim_to_mean), 
               Mean1Dim(smooth_ano_obs, dim_to_mean), 
               compROW = required_complete_row, 
               limits = c(ceiling((runmean_months + 1) / 2), 
                          leadtimes_per_startdate - floor(runmean_months / 2)))
  PlotVsLTime(corr, toptitle = "correlations", ytitle = "correlation", 
              monini = 11, limits = c(-1, 2), listexp = c('CMIP5 IC3'), 
              listobs = c('ERSST'), biglab = FALSE, hlines = c(-1, 0, 1), 
              fileout = 'tos_cor.eps')
}
\author{
  History:
  0.1  -  2011-04  (V. Guemas, \email{vguemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
