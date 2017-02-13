\name{StatSeasAtlHurr}
\alias{StatSeasAtlHurr}
\title{Compute estimate of seasonal mean of Atlantic hurricane activity}
\description{
Compute one of G. Villarini's statistically downscaled measure of mean Atlantic hurricane activity and its variance. The hurricane activity is estimated using seasonal averages of sea surface temperature anomalies over the tropical Atlantic (bounded by 10N-25N and 80W-20W) and the tropics at large (bounded by 30N-30S). The anomalies are for the JJASON season.\cr
The estimated seasonal average is either 1) number of hurricanes, 2) number of tropical cyclones with lifetime >=48h or 3) power dissipation index (PDI; in 10^11 m^3 s^{-2}).\cr
The statistical models used in this function are described in\cr
}
\usage{
StatSeasAtlHurr(atlano = NULL, tropano = NULL, hrvar = 'HR')
}
\arguments{
  \item{atlano}{
Array of Atlantic sea surface temperature anomalies. Must have the same dimension as tropano.
  }
  \item{tropano}{
Array of tropical sea surface temperature anomalies. Must have the same dimension as atlano.
  }
  \item{hrvar}{
The seasonal average to be estimated. The options are either\cr 
"HR" (hurricanes) \cr
"TC" (tropical cyclones with lifetime >=48h) \cr
"PDI" (power dissipation index) \cr
 }
}
\value{
A list composed of two matrices:\cr
  \itemize{
     \item{
1. a matrix (mean) with the seasonal average values of the desired quantity.\cr 
     }
     \item{
2. a matrix (var) of the variance of that quantity.\cr
     }
  }
The dimensions of the two matrices are the same as the dimensions of atlano/tropano.
}
\examples{
# Let AtlAno represents 5 different 5-year forecasts of seasonally averaged 
# Atlantic sea surface temperature anomalies.
AtlAno <- matrix(c(-0.31, -0.36, 0.26, -0.16, -0.16, 
                   -0.06, -0.22, -0.31, -0.36, -0.39, 
                    0.20, -0.14, 0.12, 0.22, 0.02,
                   -0.28, 0.26, -0.10, 0.18, 0.33, 
                    0.45, 0.46, 0.04, 0.12, 0.21), 
                    nrow = 5, ncol = 5)
# Let TropAno represents 5 corresponding 5-year forecasts of seasonally averaged 
# tropical sea surface temperature anomalies.
TropAno <- matrix(c(-0.22, -.13, 0.07, -0.16, -0.15,
                     0.00,  -0.03, -0.22, -0.13, -0.10,
                     0.07, -0.07, 0.17, 0.10, -0.15,
                    -0.01, 0.08, 0.07, 0.17, 0.13,
                     0.16, 0.15, -0.09, 0.03, 0.27),
                     nrow = 5, ncol = 5)
# The seasonal average of hurricanes for each of the five forecasted years, 
# for each forecast, would then be given by
hr_count <- StatSeasAtlHurr(atlano = AtlAno, 
                            tropano = TropAno, 
                            hrvar = 'HR')
print(hr_count$mean)
} 
\references{
Villarini et al. (2010) Mon Wea Rev, 138, 2681-2705.\cr
Villarini et al. (2012) Mon Wea Rev, 140, 44-65.\cr
Villarini et al. (2012) J Clim, 25, 625-637.\cr
An example of how the function can be used in hurricane forecast studies is given in\cr
Caron, L.-P. et al. (2014) Multi-year prediction skill of Atlantic hurricane activity in CMIP5 decadal hindcasts. Climate Dynamics, 42, 2675-2690. doi:10.1007/s00382-013-1773-1.
}
 \author{
History:\cr
0.1  -  2015-11  (Louis-Philippe Caron, \email{louis-philippe.caron@bsc.es})  -  Original code 
}
\keyword{datagen}
