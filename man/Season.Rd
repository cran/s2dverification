\name{Season}
\alias{Season}
\title{
  Computes Seasonal Means
}
\description{
  Computes seasonal means on timeseries organized in a matrix of any number of 
  dimensions up to 10 dimensions where the time dimension is one of those 10
  dimensions.
}
\usage{
  Season(var, posdim = 4, monini, moninf, monsup)
}
\arguments{
  \item{var}{
    Matrix containing the timeseries along one of its dimensions.
  }
  \item{posdim}{
    Dimension along which to compute seasonal means = Time dimension
  }
  \item{monini}{
    First month of the time-series: 1 to 12.
  }
  \item{moninf}{
    Month when to start the seasonal means: 1 to 12.
  }
  \item{monsup}{
    Month when to stop the seasonal means: 1 to 12.
  }
}
\value{
  Matrix with the same dimensions as var except along the posdim dimension 
  which length corresponds to the number of seasons. Partial seasons are not 
  accounted for.
}
\examples{
  startDates <- c('19901101', '19951101', '20001101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'areave', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  leadtimes_dimension <- 4
  initial_month <- 11
  mean_start_month <- 12
  mean_stop_month <- 2
  season_means_mod <- Season(sampleData$mod, leadtimes_dimension, initial_month,
                             mean_start_month, mean_stop_month)
  season_means_obs <- Season(sampleData$obs, leadtimes_dimension, initial_month,
                             mean_start_month, mean_stop_month)
  PlotAno(season_means_mod, season_means_obs, startDates, 
          toptitle = paste('winter (DJF) temperatures'), ytitle = c('K'), 
          legends = 'ERSST', biglab = FALSE, fileout = 'tos_season_means.eps')
}
\author{
  History:
  0.1  -  2011-03  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
