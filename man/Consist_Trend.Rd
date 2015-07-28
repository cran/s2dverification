\name{Consist_Trend}
\alias{Consist_Trend}
\title{
Computes Trends Using Only Model Data For Which Observations Are Available
}
\description{
Computes trends by least square fitting together with the associated error interval for both the observational and model data.\cr
Provides also the detrended observational and modeled data.\cr
The trend is computed along the second dimension, expected to be the start date dimension (the user is supposed to perform an ensemble averaging operation with \code{Mean1Dim()} prior to using \code{Consist_trend()}).
}
\usage{
Consist_Trend(var_exp, var_obs, interval = 1)
}
\arguments{
  \item{var_exp}{
Ensemble mean of model hindcasts with dimensions:\cr
  c(nmod/nexp, nsdates, nltime) up to\cr
  c(nmod/nexp, nsdates, nltime, nlevel, nlat, nlon)
  }
  \item{var_obs}{
Ensemble mean of observational data with dimensions:\cr
  c(nobs, nsdates, nltime) up to\cr
  c(nobs, nsdates, nltime, nlevel, nlat, nlon)\cr
Dimensions 2 to 6 should be the same as var_exp.
  }
  \item{interval}{
Number of months/years between 2 start dates. Default = 1. The trends will be provided respectively in field unit per month or per year.
  }
}
\value{
  \item{$trend}{
Trends of model and observational data with dimensions:\cr
  c(nmod/nexp + nobs, 3, nltime) up to\cr
  c(nmod/nexp + nobs, 3, nltime, nlevel, nlat, nlon)\cr
The length 3 dimension corresponds to the lower limit of the 95\% confidence interval, the slope of the trends and the upper limit of the 95\% confidence interval.
  }
  \item{$detrendedmod}{
Same dimensions as var_exp with linearly detrended values of var_exp along the second = start date dimension.
  }
  \item{$detrendedobs}{
Same dimensions as var_exp with linearly detrended values of var_obs along the second = start date dimension.
  }
}
\examples{
# Load sample data as in Load() example:
example(Load)
clim <- Clim(sampleData$mod, sampleData$obs)
ano_exp <- Ano(sampleData$mod, clim$clim_exp)
ano_obs <- Ano(sampleData$obs, clim$clim_obs)
runmean_months <- 12
dim_to_smooth <- 4  # Smooth along lead-times
smooth_ano_exp <- Smoothing(ano_exp, runmean_months, dim_to_smooth)
smooth_ano_obs <- Smoothing(ano_obs, runmean_months, dim_to_smooth)
dim_to_mean <- 2  # Mean along members
years_between_startdates <- 5
trend <- Consist_Trend(Mean1Dim(smooth_ano_exp, dim_to_mean), 
                       Mean1Dim(smooth_ano_obs, dim_to_mean), 
                       years_between_startdates)

PlotVsLTime(trend$trend, toptitle = "trend", ytitle = "K/(5 years)", 
            monini = 11, limits = c(-0.8, 0.8), listexp = c('CMIP5 IC3'), 
            listobs = c('ERSST'), biglab = FALSE, hlines = c(0), 
            fileout = 'tos_consist_trend.eps')
PlotAno(InsertDim(trend$detrendedmod,2,1), InsertDim(trend$detrendedobs,2,1), 
        startDates, "Detrended tos anomalies", ytitle = 'K', 
        legends = 'ERSST', biglab = FALSE, fileout = 'tos_detrended_ano.eps')
}
\author{
History:\cr
0.1  -  2011-11  (V. Guemas, \email{vguemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
