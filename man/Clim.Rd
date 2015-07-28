\name{Clim}
\alias{Clim}
\title{Computes Per-pair/Kharin/Fuckar Climatologies}
\description{
This function computes climatologies from the experimental and observational matrices output from \code{Load()} using one of the following methods:
  \itemize{
    \item{
1) per-pair method (Garcia-Serrano and Doblas-Reyes, CD, 2012)
    }
    \item{
2) Kharin method (Karin et al, GRL, 2012)
    }
    \item{
3) Fuckar method (Fuckar et al, GRL, 2014)
    }
  }
Clim function computes climatologies using the startdates covered by the whole experiments/observational data sets. The startdates not available for all the data (model and obs) are excluded when computing the climatologies.
}
\usage{
Clim(var_exp, var_obs, memb = TRUE, kharin = FALSE, NDV = FALSE)
}
\arguments{
  \item{var_exp}{
Model data:\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime) up to\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)
  }
  \item{var_obs}{
Observational data: \cr
  c(nobs, nmemb, nsdates, nltime) up to\cr
  c(nobs, nmemb, nsdates, nltime, nlevel, nlat, nlon)
  }
  \item{memb}{
memb: TRUE/FALSE (1 climatology for each member). Default = TRUE.
  }
  \item{kharin}{
TRUE/FALSE (if Kharin method is applied or not). Default = FALSE.
  }
  \item{NDV}{
TRUE/FALSE (if Fuckar method is applied or not). Default = FALSE. 
  }
}
\value{
  \item{clim_exp}{Array with same dimensions as var_exp except the third (starting dates) and, depending on the parameters, the second (members), which disappear.}
  \item{clim_obs}{Array with same dimensions as var_obs except the third (starting dates) and, depending on the parameters, the second (members), which disappear.}
}
\examples{
# Load sample data as in Load() example:
example(Load)
clim <- Clim(sampleData$mod, sampleData$obs)
PlotClim(clim$clim_exp, clim$clim_obs, 
         toptitle = paste('sea surface temperature climatologies'), 
         ytitle = 'K', monini = 11, listexp = c('CMIP5 IC3'), 
         listobs = c('ERSST'), biglab = FALSE, fileout = 'tos_clim.eps')
}
\author{
History:\cr
0.9  -  2011-03  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
