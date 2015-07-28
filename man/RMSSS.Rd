\name{RMSSS}
\alias{RMSSS}
\title{
Computes Root Mean Square Skill Score
}
\description{
Arrays var_exp & var_obs should have the same dimensions except along posloop where the length can be different, with the number of experiments/models for var_exp (nexp) and the number of obserational datasets for var_obs (nobs).\cr
RMSSS computes the Root Mean Square Skill Score of each jexp in 1:nexp against each jobs in 1:nobs which gives nexp x nobs RMSSS for each other grid point of the matrix (each latitude/longitude/level/leadtime).\cr
The RMSSS are computed along the posRMS dimension which should correspond to the startdate dimension.\cr
The p-value is provided by a one-sided Fisher test.
}
\usage{
RMSSS(var_exp, var_obs, posloop = 1, posRMS = 2)
}
\arguments{
  \item{var_exp}{
Array of experimental data.
  }
  \item{var_obs}{
Array of observational data, same dimensions as var_exp except along posloop dimension, where the length can be nobs instead of nexp.
  }
  \item{posloop}{
Dimension nobs and nexp.
  }
  \item{posRMS}{
Dimension along which the RMSE are to be computed (the dimension of the start dates).
  }
}
\value{
Array with dimensions :\cr
  c(length(posloop) in var_exp, length(posloop) in var_obs, 2, all other dimensions of var_exp & var_obs except posRMS).\cr
The dimension 2 corresponds to the RMSSS and the p.value of the one-sided 
Fisher test with Ho: RMSSS = 0.
}
\examples{
# Load sample data as in Load() example:
example(Load)
clim <- Clim(sampleData$mod, sampleData$obs)
ano_exp <- Ano(sampleData$mod, clim$clim_exp)
ano_obs <- Ano(sampleData$obs, clim$clim_obs)
rmsss <- RMSSS(Mean1Dim(ano_exp, 2), Mean1Dim(ano_obs, 2))
rmsss2 <- array(dim = c(dim(rmsss)[1:2], 4, dim(rmsss)[4]))
rmsss2[, , 2, ] <- rmsss[, , 1, ]
rmsss2[, , 4, ] <- rmsss[, , 2, ]
PlotVsLTime(rmsss, toptitle = "Root Mean Square Skill Score", ytitle = "", 
            monini = 11, limits = c(-1, 1.3), listexp = c('CMIP5 IC3'), 
            listobs = c('ERSST'), biglab = FALSE, hlines = c(-1, 0, 1), 
            fileout = 'tos_rmsss.eps')
}
\author{
History:\cr
0.1  -  2012-04  (V. Guemas, \email{vguemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
