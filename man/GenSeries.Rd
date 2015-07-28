\name{GenSeries}
\alias{GenSeries}
\title{Generates An AR1 Time Series}
\description{This functions generates AR1 processes containing n data, with alpha as autocorrelation at lag 1, and mean and standard deviation provided by the mean and std arguments.}
\usage{
GenSeries(n, alpha, mean, std)
}
\arguments{
  \item{n}{Length of the timeseries to be generated.}
  \item{alpha}{Autocorrelation at lag 1.}
  \item{mean}{Mean of the data.}
  \item{std}{Standard deviation of the data.}
}
\value{AR1 timeseries}
\examples{
series <- GenSeries(1000, 0.35, 2, 1)
plot(series, type = 'l')
}
\author{
History:\cr
0.1  -  2012-04  (L. Auger, \email{ludovic.auger@meteo.fr})  -  Original code\cr
1.0  -  2012-04  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
