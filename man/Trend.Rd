\name{Trend}
\alias{Trend}
\title{
Computes Trends
}
\description{
Computes the trend along the posTR dimension of the matrix var by least square fitting, and the associated an error interval.\cr
Provide also the detrended data.\cr
The confidence interval relies on a student-T distribution.
}
\usage{
Trend(var, posTR = 2, interval = 1)
}
\arguments{
  \item{var}{
Matrix of any number of dimensions up to 10.
  }
  \item{posTR}{
Position along which to compute the trend.
  }
  \item{interval}{
Number of months/years between 2 points along posTR dimension.\cr
Default = 1.\cr
The trend would be provided in number of units per month or year.
  }
}
\value{
  \item{$trend}{
Same dimensions as var except along the posTR dimension which is replaced by a length 3 dimension, corresponding to the lower limit of the  95\% confidence interval, trends and the upper limit of the 95\% confidence interval for each point of the matrix along all the other dimensions.
  }
  \item{$detrended}{
Same dimensions as var with linearly detrended var along the posTR dimension.
  }
}
\examples{
# Load sample data as in Load() example:
example(Load)
months_between_startdates <- 60
trend <- Trend(sampleData$obs, 3, months_between_startdates)
PlotVsLTime(trend$trend, toptitle = "trend", ytitle = "K / (5 year)",
            monini = 11, limits = c(-1,1), listexp = c('CMIP5 IC3'),
            listobs = c('ERSST'), biglab = FALSE, hlines = 0,
            fileout = 'tos_obs_trend.eps')
PlotAno(trend$detrended, NULL, startDates, 
        toptitle = 'detrended anomalies (along the startdates)', ytitle = 'K',
        legends = 'ERSST', biglab = FALSE, fileout = 'tos_detrended_obs.eps')
}
\author{
History:\cr
0.1  -  2011-05  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
