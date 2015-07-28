\name{Ano}
\alias{Ano}
\title{Computes Forecast or Observed Anomalies}
\description{
This function computes anomalies from any experimental or observational matrix output from \code{Load()} and their climatologies output from \code{Clim()}.
}
\usage{Ano(var, clim)}
\arguments{
  \item{var}{
Model or observational data:\cr
  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime) up to\cr
  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)
  }
  \item{clim}{
Climatologies from clim:\cr
  c(nmod/nexp/nobs, nltime) up to\cr
  c(nmod/nexp/nobs, nltime, nlevel, nlat, nlon)\cr
or\cr
  c(nmod/nexp/nobs, nmemb/nparam, nltime) up to\cr
  c(nmod/nexp/nobs, nmemb/nparam, nltime, nlevel, nlat, nlon)\cr
or\cr
  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime) up to\cr
  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)\cr
depending on the options provided to \code{Clim()}
  }
}
\value{
Array with same dimensions as 'var'.
}
\examples{
# Load sample data as in Load() example:
example(Load)
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
History:\cr
0.1  -  2012-03  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
