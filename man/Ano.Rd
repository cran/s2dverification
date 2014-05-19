\name{Ano}
\alias{Ano}
\title{Computes Forecast or Observed Anomalies}
\description{
  This function computes anomalies from any experimental or observational 
  matrix output from \code{Load()} and their climatologies output from 
  \code{Clim()}.
}
\usage{
  Ano(var, clim)
}
\arguments{
  \item{var}{
    Model or observational data:
      c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime) up to
      c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)
  }
  \item{clim}{
    Climatologies from clim:
      c(nmod/nexp/nobs, nltime) up to
      c(nmod/nexp/nobs, nltime, nlevel, nlat, nlon)
      or
      c(nmod/nexp/nobs, nmemb/nparam, nltime) up to
      c(nmod/nexp/nobs, nmemb/nparam, nltime, nlevel, nlat, nlon)
      or
      c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime) up to
      c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)
      depending on the options provided to \code{Clim()}
  }
}
\value{
  Matrix with same dimensions as var
}
\examples{
  startDates <- c('19901101', '19951101', '20001101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'areave', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  clim <- Clim(sampleData$mod, sampleData$obs)
  ano_exp <- Ano(sampleData$mod, clim$clim_exp)
  ano_obs <- Ano(sampleData$obs, clim$clim_obs)
  runmean_nb_months <- 12
  dim_to_smooth <- 4  # Smooth along lead-times
  smooth_ano_exp <- Smoothing(ano_exp, runmean_nb_months, dim_to_smooth)
  smooth_ano_obs <- Smoothing(ano_obs, runmean_nb_months, dim_to_smooth)
  PlotAno(smooth_ano_exp, smooth_ano_obs, startDates, 
          toptitle = paste('smoothed anomalies'), ytitle = c('K', 'K', 'K'), 
          legends = 'ERSST', biglab = FALSE, fileout = 'tos_ano.eps')
}
\author{
  History:
  0.1  -  2012-03  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
