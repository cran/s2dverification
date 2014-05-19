\name{ACC}
\alias{ACC}
\title{
  Computes Anomaly Correlation Coefficient (Spatial Correlation)
}
\description{
  Matrices var_exp and var_obs should have dimensions:
    c(nexp/nobs, nsdates, nltimes, nlat, nlon).
  The dimension 2 to 6 should have the same lengths.
  ACC computes the Anomaly Correlation Coefficient for each jexp in 1:nexp and 
  each jobs in 1:nobs which gives nexp x nobs ACC for each startdate and each 
  leadtime.
  The confidence interval is computed by a Fisher transformation.
  The significance level relies on a one-sided student-T distribution.
  A domain can be selected by providing the longitudes/latitudes (lon/lat) of 
  the grid together with the corner of the domain to be selected:
    lonlatbox = c(lonmin, lonmax, latmin, latmax)
}
\usage{
  ACC(var_exp, var_obs, lon = NULL, lat = NULL, lonlatbox = NULL)
}
\arguments{
  \item{var_exp}{
    Matrix of experimental data with dimensions: 
      c(nexp, nsdates, nltimes, nlat, nlon)
  }
  \item{var_obs}{
    Matrix of observational data with the same dimensions as var_exp except 
    along the first dimension:
      c(nobs, nsdates, nltimes, nlat, nlon)
  }
  \item{lon}{
    Array of longitudes of the var_exp/var_obs grids, optional.
  }
  \item{lat}{
    Array of latitudes of the var_exp/var_obs grids, optional.
  }
  \item{lonlatbox}{
    Domain to select : c(lonmin, lonmax, latmin, latmax), optional.
  }
}
\value{
  \item{ACC}{
    Matrix with dimensions : c(nexp, nobs, nsdates, nleadtimes, 4) 
    The fifth dimension of length 4 corresponds to the lower limit of the 95\% 
    confidence interval, the ACC, the upper limit of the 95\% confidence 
    interval and the 95\% significance level. 
  }
  \item{MACC}{
    Mean Anomaly Correlation Coefficient with dimensions :
      c(nexp, nobs, nleadtimes) 
  }
}
\examples{
  startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 4, output = 'lonlat', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  sampleData$mod <- Season(sampleData$mod, 4, 11, 12, 2)
  sampleData$obs <- Season(sampleData$obs, 4, 11, 12, 2)
  clim <- Clim(sampleData$mod, sampleData$obs)
  ano_exp <- Ano(sampleData$mod, clim$clim_exp)
  ano_obs <- Ano(sampleData$obs, clim$clim_obs)
  acc <- ACC(Mean1Dim(sampleData$mod, 2), 
             Mean1Dim(sampleData$obs, 2))
  PlotACC(acc$ACC, startDates)
}
\author{
  History:
  0.1  -  2013-08  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
