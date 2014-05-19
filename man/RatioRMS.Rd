\name{RatioRMS}
\alias{RatioRMS}
\title{
  Computes The Ratio Between The RMSE Scores of 2 Experiments.
}
\description{
  Matrix var_exp1 / var_exp2 / var_obs should have the same dimensions.
  The ratio RMSE(var_exp1, var_obs) / RMSE(var_exp2, var_obs) is output.
  The p-value is provided by a two-sided Fischer test.
}
\usage{
  RatioRMS(var_exp1, var_exp2, var_obs, posRMS = 1)
}
\arguments{
  \item{var_exp1}{
    Matrix of experimental data 1.
  }
  \item{var_exp2}{
    Matrix of experimental data 2, same dimensions as var_exp1.
  }
  \item{var_obs}{
    Matrix of observational data, same dimensions as var_exp1.
  }
  \item{posRMS}{
    Dimension along which the RMSE are to be computed = the position of the
    start dates.
  }
}
\value{
  Matrix with the same dimensions as var_exp1/var_exp2/var_obs except along 
  posRMS where the dimension has length 2.
  The dimension 2 corresponds to the ratio between the RMSE (RMSE1/RMSE2) and 
  the p.value of the two-sided Fisher test with Ho: RMSE1/RMSE2 = 1.
}
\examples{
  startDates <- c('19901101', '19951101', '20001101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'lonlat', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  leadtimes_dimension <- 4
  initial_month <- 11
  mean_start_month <- 12
  mean_stop_month <- 2                                
  sampleData$mod <- Season(sampleData$mod, leadtimes_dimension, initial_month, 
                           mean_start_month, mean_stop_month)
  sampleData$obs <- Season(sampleData$obs, leadtimes_dimension, initial_month, 
                           mean_start_month, mean_stop_month)                              
  clim <- Clim(sampleData$mod, sampleData$obs)
  ano_exp <- Ano(sampleData$mod, clim$clim_exp)
  ano_obs <- Ano(sampleData$obs, clim$clim_obs)
  rrms <- RatioRMS(Mean1Dim(ano_exp[ , 1:3, , , , ], 1)[, 1, , ], 
                   Mean1Dim(ano_exp[ , 4:5, , , , ], 1)[, 1, , ],
                   Mean1Dim(ano_obs, 2)[1, , 1, , ], 1)
  PlotEquiMap(rrms[1, , ], sampleData$lon, sampleData$lat, 
              toptitle = 'Ratio RMSE')
}
\author{
  History:
  0.1  -  2011-11  (V. Guemas, \email{vguemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
