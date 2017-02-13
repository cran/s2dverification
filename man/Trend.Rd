\name{Trend}
\alias{Trend}
\alias{.Trend}
\title{
Computes the Trend of the Ensemble Mean
}
\description{
Computes the trend along the forecast time of the ensemble mean by least square fitting, and the associated error interval.\cr
Trend() also provides the time series of the detrended ensemble mean forecasts.\cr
The confidence interval relies on a student-T distribution.

.Trend provides the same functionality but taking a matrix ensemble members as input (exp).
}
\usage{
Trend(var, posTR = 2, interval = 1, siglev = 0.95, conf = TRUE)

.Trend(exp, interval = 1, siglev = 0.95, conf = TRUE)
}
\arguments{
 \item{var}{
Array of any number of dimensions up to 10.
  }
  \item{exp}{
M by N matrix of M forecasts from N ensemble members.
  } 
  \item{interval}{
Number of months/years between 2 points along posTR dimension.\cr
Default = 1.\cr
The trend would be provided in number of units per month or year.
  }
  \item{siglev}{
Confidence level for the computation of confidence interval. 0.95 by default.
  }
  \item{conf}{
Whether to compute the confidence levels or not. TRUE by default.
  }
   \item{posTR}{
Position along which to compute the trend.
  }
}
\value{
  \item{$trend}{
The intercept and slope coefficients for the least squares fitting of the trend, 
  }
  \item{$conf.int}{
Corresponding to the limits of the \code{siglev}\% confidence interval (only present if \code{conf = TRUE}) for the slope coefficient.
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
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN\cr
2.0  -  2017-02  (A. Hunter, \email{alasdair.hunter at bsc.es})  -  Adapt to veriApply()
}
\keyword{datagen}
