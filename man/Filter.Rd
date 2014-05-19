\name{Filter}
\alias{Filter}
\title{Filter Frequency Peaks From An Array}
\description{
This function filters from the xdata array, the signal of frequency freq.\cr 
The filtering is performed by dichotomy, seeking for the frequency around freq and the phase that maximizes the signal to subtract to xdata.\cr
The maximization of the signal to subtract relies on a minimization of the mean square differences between xdata and a cosine of given frequency and phase.
}
\usage{
Filter(xdata, freq)
}
\arguments{
  \item{xdata}{Array to be filtered.}
  \item{freq}{Frequency to filter.}
}
\value{Filtered Array}
\examples{
  startDates <- c('19901101', '19951101', '20001101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'areave', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)

  for (jmemb in 1:5) {
    for (jstartdate in 1:3) {
      spectrum <- Spectrum(sampleData$mod[1, jmemb, jstartdate, ])
      for (jlen in 1:dim(spectrum)[1]) {
        if (spectrum[jlen, 2] > spectrum[jlen, 4]) {
          sampleData$mod[1, jmemb, jstartdate, ] <- Filter(sampleData$mod[1, jmemb, jstartdate, ], 
                                                           spectrum[jlen, 1])
        }
      }
    }
  }

  PlotAno(sampleData$mod, sdates = startDates, fileout = 'filtered_data.eps') 
}
\author{
History:\cr
0.1  -  2012-02  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
1.0  -  2012-02  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
