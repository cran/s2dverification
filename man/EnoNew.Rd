\name{EnoNew}
\alias{EnoNew}
\title{Computes Effective Sample Size Following Guemas et al, BAMS, 2013b}
\description{This function computes the equivalent number of independent data 
  in the xdata array following the method described in Guemas V., Auger L., 
  Doblas-Reyes F., JAMC, 2013. The method relies on the Trenberth (1984) 
  formula combined with a reduced uncertainty of the estimated autocorrelation 
  function compared to the original approach.}
\usage{
EnoNew(xdata, detrend = FALSE, filter = FALSE)
}
\arguments{
  \item{xdata}{Timeseries from which the equivalent number of independent data is requested}
  \item{detrend}{TRUE applies a linear detrending to xdata prior to the estimation of the equivalent number of independant data.}
  \item{filter}{TRUE applies a filtering of any frequency peak prior to the estimation of the equivalent number of independant data.}
}
\examples{
  startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'lonlat', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  eno <- EnoNew(sampleData$mod[1, 1, , 1, 10, 15])
  print(eno)
}
\author{
History:\cr
0.1  -  2012-06  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
