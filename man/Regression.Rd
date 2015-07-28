\name{Regression}
\alias{Regression}
\title{
Computes The Regression Of A Matrix On Another Along A Dimension
}
\description{
Computes the regression of the input matrice vary on the input matrice varx along the posREG dimension by least square fitting. Provides the slope of the regression, the associated confidence interval, and the intercept.\cr
Provides also the vary data filtered out from the regression onto varx.\cr
The confidence interval relies on a student-T distribution.
}
\usage{
Regression(vary, varx, posREG = 2)
}
\arguments{
  \item{vary}{
Matrix of any number of dimensions up to 10.
  }
  \item{varx}{
Matrix of any number of dimensions up to 10. Same dimensions as vary.
  }
  \item{posREG}{
Position along which to compute the regression.
  }
}
\value{
  \item{$regression}{
Matrix with same dimensions as varx and vary except along posREG dimension which is replaced by a length 4 dimension, corresponding to the lower limit of the 95\% confidence interval, the slope, the upper limit of the 95\% confidence interval and the intercept.
  }
  \item{$filtered}{
Same dimensions as vary filtered out from the regression onto varx along the posREG dimension.
  }
}
\examples{
# See examples on Load() to understand the first lines in this example
  \dontrun{
configfile <- paste0(tempdir(), '/sample.conf')
ConfigFileCreate(configfile, confirm = FALSE)
c <- ConfigFileOpen(configfile)
c <- ConfigEditDefinition(c, 'DEFAULT_GRID', 'r20x10', confirm = FALSE)
c <- ConfigEditDefinition(c, 'DEFAULT_VAR_MIN', '-1e19', confirm = FALSE)
c <- ConfigEditDefinition(c, 'DEFAULT_VAR_MAX', '1e19', confirm = FALSE)
c <- ConfigAddVar(c, '2d', 'tos')
data_path <- system.file('sample_data', package = 's2dverification')
exp_data_path <- paste0(data_path, '/model/$EXP_NAME$')
obs_data_path <- paste0(data_path, '/$OBS_NAME$')
c <- ConfigAddEntry(c, 'experiments', 'file-per-startdate',
     dataset_name = 'experiment', var_name = 'tos', main_path = exp_data_path,
     file_path = '$STORE_FREQ$_mean/$VAR_NAME$_3hourly/$VAR_NAME$_$START_DATE$.nc')
c <- ConfigAddEntry(c, 'observations', 'file-per-month',
     dataset_name = 'observation', var_name = 'tos', main_path = obs_data_path,
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
sampleData$mod <- Season(sampleData$mod, 4, 11, 12, 2)
sampleData$obs <- Season(sampleData$obs, 4, 11, 12, 2)
reg <- Regression(Mean1Dim(sampleData$mod, 2),
                  Mean1Dim(sampleData$obs, 2), 2)
PlotEquiMap(reg$regression[1, 2, 1, , ], sampleData$lon, sampleData$lat, 
            toptitle='Regression of the prediction on the observations', 
            sizetit = 0.5)
}
\author{
History:\cr
0.1  -  2013-05  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
