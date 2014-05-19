\name{Histo2Hindcast}
\alias{Histo2Hindcast}
\title{
  Chunks Long Simulations For Comparison With Hindcasts
}
\description{
  This function reorganizes a long run (historical typically) with only one 
  start date into chunks corresponding to a set of start dates. The expected 
  input structure is the one output from \code{Load()} with 4 to 7 dimensions.
}
\usage{
  Histo2Hindcast(varin, sdatesin, sdatesout, nleadtimesout)
}
\arguments{
  \item{varin}{
    Input model or observational data:
      c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltimes) up to
      c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltimes, nlevel, nlat, nlon)
  }
  \item{sdatesin}{
    Start date of the input matrix 'YYYYMMDD'.
  }
  \item{sdatesout}{
    List of start dates of the output matrix c('YYYYMMDD', 'YYYYMMDD', ...).
  }
  \item{nleadtimesout}{
    Number of leadtimes in the output matrix.
  }
}
\value{
  A matrix with the same number of dimensions as the input one, the same 
  dimensions 1 and 2 and potentially the same dimensions 5 to 7. Dimensions 3 
  and 4 are set by the arguments sdatesout and nleadtimesout.
}
\examples{
  startDates <- c('19901101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'areave', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  start_dates_out <- c('19901101', '19911101', '19921101', '19931101', '19941101')
  leadtimes_per_startdate <- 12
  experimental_data <- Histo2Hindcast(sampleData$mod, startDates[1], 
                                      start_dates_out, leadtimes_per_startdate)
  observational_data <- Histo2Hindcast(sampleData$obs, startDates[1], 
                                       start_dates_out, leadtimes_per_startdate)
  PlotAno(experimental_data, observational_data, start_dates_out, 
          toptitle = paste('anomalies reorganized into shorter chunks'), 
          ytitle = 'K', fileout='tos_histo2hindcast.eps')
}
\author{
  History:
  0.1  -  2012-11  (V. Guemas, \email{vguemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
