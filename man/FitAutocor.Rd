\name{FitAutocor}
\alias{FitAutocor}
\title{Fits an AR1 Autocorrelation Function Using Dichotomy}
\description{This function fits the theoretical autocorrelation function of an AR1 to the first guess of estimated autocorrelation function estacf containing any number of lags. The fitting relies on a dichotomial minimisation of the mean square differences between both autocorrelation functions. It returns the autocorrelation at lag 1 of the fitted AR1 process.}
\usage{
FitAutocor(estacf, window = c(-1, 1), prec = 0.01)
}
\arguments{
  \item{estacf}{First guess of the autocorrelation function}
  \item{window}{Interval in which the autocorrelation at lag 1 should be found.}
  \item{prec}{Precision to which the autocorrelation function at lag 1 is to be estimated.}
}
\value{Best estimate of the autocorrelation at lag 1}
\examples{
series <- GenSeries(1000, 0.35, 2, 1)
estacf <- acf(series[951:1000], plot = FALSE)$acf
alpha <-  FitAutocor(estacf, c(-1, 1), 0.01)
print(alpha)
}
\author{
History:\cr
0.1  -  2012-02  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
