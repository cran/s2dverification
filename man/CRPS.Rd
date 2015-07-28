\name{CRPS}
\alias{CRPS}
\title{
Compute Continuous Ranked Probability Score (CRPS) For Ensemble Forecasts
}
\description{
Returns the value of the CRPS.
}
\usage{
CRPS(obs, pred)
}
\arguments{
  \item{obs}{
Vector of observations for a continuous variable (e.g. rainfall)
  }
  \item{pred}{
Matrix of ensemble forecasts for a continuous variable (e.g. rainfall).\cr
It should have the following structure pred[nmembers,nsdates]
  }
  
}
\value{
  \item{CRPS}{
Continuous ranked probability score
    }
  }
\references{
Wilks (2006) Statistical Methods in the Atmospheric Sciences.\cr
Jolliffe and Stephenson (2012) Forecast verification: A practitioner's guide in the Atmospheric science\cr
Hersbach (2000) Weather and Forecasting 15:559-570
}
\examples{
a <- runif(10)
b <- array(, dim = c(10,10))
for (i in 1:10) {
  b[i, ] <- runif(10)
}
x <- CRPS(a, b)
x$CRPS
}
\author{
History:\cr
1.0  -  2013-03-12  (L. Rodrigues, \email{lrodrigues at ic3.cat})  -  Original code\cr
1.1  -  2014-11-21(V. Torralba,\email{veronica.torralba at ic3.cat})  -  Changes in the documentation
}
\keyword{datagen}
