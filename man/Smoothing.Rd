\name{Smoothing}
\alias{Smoothing}
\title{
  Smoothes A Matrix Along A Dimension
}
\description{
  Smoothes a matrix of any number of dimensions up to 10 dimensions along one of its dimensions
}
\usage{
  Smoothing(var, runmeanlen = 12, numdimt = 4)
}
\arguments{
  \item{var}{
    Matrix to be smoothed along one of its dimension (typically the forecast time dimension).
  }
  \item{runmeanlen}{
    Running mean length in number of sampling units (typically months).
  }
  \item{numdimt}{
    Dimension to smooth.
  }
}
\value{
  Matrix with same the dimensions as var but smoothed along the numdimt
  dimension.
}
\examples{
  startDates <- c('19901101', '19951101', '20001101')
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
  PlotAno(smooth_ano_exp, smooth_ano_obs, startDates, 
          toptitle = "Smoothed Mediterranean mean SST", ytitle = "K",
          fileout = "tos_smoothed_ano.eps")
}
\author{
  History:
  0.1  -  2011-03  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
