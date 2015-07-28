\name{RatioSDRMS}
\alias{RatioSDRMS}
\title{
Computes The Ratio Between the Ensemble Spread and the RMSE of the Ensemble Mean
}
\description{
Matrices var_exp & var_obs should have dimensions between\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime)\cr
and\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)\cr
The ratio between the standard deviation of the members around the ensemble mean in var_exp and the RMSE between var_exp and var_obs is output for each experiment and each observational dataset.\cr
The p-value is provided by a one-sided Fischer test.
}
\usage{
RatioSDRMS(var_exp, var_obs)
}
\arguments{
  \item{var_exp}{
Model data:\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime) up to\cr
  c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)
  }
  \item{var_obs}{
Observational data:\cr
  c(nobs, nmemb, nsdates, nltime) up to\cr
  c(nobs, nmemb, nsdates, nltime, nlevel, nlat, nlon)
  }
}
\value{
Matrix with dimensions c(nexp/nmod, nobs, 2, nltime) up to
  c(nexp/nmod, nobs, 2, nltime, nlevel, nlat, nlon) dimensions.
The dimension 2 corresponds to the ratio (SD/RMSE) and the p.value of the one-sided Fisher test with Ho: SD/RMSE = 1.
}
\examples{
# Load sample data as in Load() example:
example(Load)
rsdrms <- RatioSDRMS(sampleData$mod, sampleData$obs)
rsdrms2 <- array(dim = c(dim(rsdrms)[1:2], 4, dim(rsdrms)[4]))
rsdrms2[, , 2, ] <- rsdrms[, , 1, ]
rsdrms2[, , 4, ] <- rsdrms[, , 2, ]
PlotVsLTime(rsdrms2, toptitle = "Ratio ensemble spread / RMSE", ytitle = "", 
            monini = 11, limits = c(-1, 1.3), listexp = c('CMIP5 IC3'), 
            listobs = c('ERSST'), biglab = FALSE, siglev = TRUE, 
            fileout = 'tos_rsdrms.eps')
}
\author{
History:\cr
0.1  -  2011-12  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau-manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
