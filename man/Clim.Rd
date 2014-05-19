\name{Clim}
\alias{Clim}
\title{Computes Per-pair/Kharin/FALSEuckar Climatologies}
\description{
  This function computes climatologies from the experimental and 
  observational matrices output from \code{Load()} using one of the 
  following methods:
  1) per-pair method (Garcia-Serrano and Doblas-Reyes, CD, 2012)
  2) Kharin method (Karin et al, GRL, 2012)
  3) Fuckar mathod (Fuckar et al, GRL, 2014)
}
\usage{
  Clim(var_exp, var_obs, memb = TRUE, kharin = FALSE, NDV = FALSE)
}
\arguments{
  \item{var_exp}{
    Model data:
      c(nmod/nexp, nmemb/nparam, nsdates, nltime) up to
      c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)
  }
  \item{var_obs}{
    Observational data: 
      c(nobs, nmemb, nsdates, nltime) up to
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
  \item{clim_exp}{Matrix with same dimensions as var_exp}
  \item{clim_obs}{Matrix with same dimensions as var_obs}
}
\examples{
  startDates <- c('19901101', '19951101', '20001101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'areave', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  clim <- Clim(sampleData$mod, sampleData$obs)
  PlotClim(clim$clim_exp, clim$clim_obs, 
           toptitle = paste('sea surface temperature climatologies'), 
           ytitle = 'K', monini = 11, listexp = c('CMIP5 IC3'), 
           listobs = c('ERSST'), biglab = FALSE, fileout = 'tos_clim.eps')
}
\author{
  History:
  0.9  -  2011-03  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
