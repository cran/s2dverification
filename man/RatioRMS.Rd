\name{RatioRMS}
\alias{RatioRMS}
\title{
Computes The Ratio Between The RMSE Scores of 2 Experiments.
}
\description{
Matrix var_exp1 / var_exp2 / var_obs should have the same dimensions.\cr
The ratio RMSE(var_exp1, var_obs) / RMSE(var_exp2, var_obs) is output.\cr
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
Dimension along which the RMSE are to be computed = the position of the start dates.
  }
}
\value{
Matrix with the same dimensions as var_exp1/var_exp2/var_obs except along posRMS where the dimension has length 2.\cr
The dimension 2 corresponds to the ratio between the RMSE (RMSE1/RMSE2) and the p.value of the two-sided Fisher test with Ho: RMSE1/RMSE2 = 1.
}
\examples{
# See examples on Load() to understand the first lines in this example
  \dontrun{
configfile <- paste0(tempdir(), '/sample.conf')
ConfigFileCreate(configfile, confirm = FALSE)
c <- ConfigFileOpen(configfile)
c <- ConfigEditDefinition(c, 'DEFAULT_VAR_MIN', '-1e19', confirm = FALSE)
c <- ConfigEditDefinition(c, 'DEFAULT_VAR_MAX', '1e19', confirm = FALSE)
data_path <- system.file('sample_data', package = 's2dverification')
exp_data_path <- paste0(data_path, '/model/$EXP_NAME$/')
obs_data_path <- paste0(data_path, '/$OBS_NAME$/')
c <- ConfigAddEntry(c, 'experiments', dataset_name = 'experiment', 
     var_name = 'tos', main_path = exp_data_path,
     file_path = '$STORE_FREQ$_mean/$VAR_NAME$_3hourly/$VAR_NAME$_$START_DATE$.nc')
c <- ConfigAddEntry(c, 'observations', dataset_name = 'observation', 
     var_name = 'tos', main_path = obs_data_path,
     file_path = '$STORE_FREQ$_mean/$VAR_NAME$/$VAR_NAME$_$YEAR$$MONTH$.nc')
ConfigFileSave(c, configfile, confirm = FALSE)

# Now we are ready to use Load().
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- Load('tos', c('experiment'), c('observation'), startDates, 
                   output = 'lonlat', latmin = 27, latmax = 48, lonmin = -12,
                   lonmax = 40, configfile = configfile)
  }
  \dontshow{
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'), 
                                                c('observation'), startDates,
                                                output = 'lonlat', 
                                                latmin = 27, latmax = 48, 
                                                lonmin = -12, lonmax = 40) 
  }
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
rrms <- RatioRMS(Mean1Dim(ano_exp[ , 1:2, , , , ], 1)[, 1, , ], 
                 ano_exp[ , 3, , , , ][, 1, , ],
                 Mean1Dim(ano_obs, 2)[1, , 1, , ], 1)
PlotEquiMap(rrms[1, , ], sampleData$lon, sampleData$lat, 
            toptitle = 'Ratio RMSE')
}
\author{
History:\cr
0.1  -  2011-11  (V. Guemas, \email{vguemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
