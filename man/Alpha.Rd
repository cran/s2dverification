\name{Alpha}
\alias{Alpha}
\title{Estimates AutoCorrelation At Lag 1 following Guemas et al, BAMS, 2013b}
\description{This function, relying on the \code{FitAcfCoef()} function, estimates the autocorrelation at lag 1 of the xdata array following the method described  in Guemas V., Auger L., Doblas-Reyes F., JAMC, 2013. After applying a linear
detrending and/or a filtering of any frequency peak if requested, the sample autocorrelation is estimated. Then the theoretical autocorrelation of an AR1 is fitted to the sample autocorrelation using the Cardano's formula (see \code{FitAcfCoef()}) to obtain the autocorrelation at lag 1. This method assumes xdata is an AR1 process.}
\usage{Alpha(xdata, detrend = FALSE, filter = FALSE)}
\arguments{
  \item{xdata}{Timeseries from which the autocorrelation at lag 1 is requested.}
  \item{detrend}{TRUE applies a linear detrending to xdata prior to the estimation of the autocorrelation at lag 1.}
  \item{filter}{TRUE applies a filtering of any frequency peak prior to the estimation of the autocorrelation at lag 1.}
}
\value{Autocorrelation at lag 1}
\examples{
  startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'areave', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  alpha <- Alpha(sampleData$mod[1, 1, , 1])
  print(alpha)
}
\author{
History:\cr
0.1  -  2012-06  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
