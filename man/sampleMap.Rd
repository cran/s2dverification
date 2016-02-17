\name{sampleMap}
\docType{data}
\alias{sampleMap}
\title{Sample Of Observational And Experimental Data For Forecast Verification In Function Of Longitudes And Latitudes}
\description{
This data set provides data in function of longitudes and latitudes for the variable 'tos', i.e. sea surface temperature, over the mediterranean zone from the sample experimental and observational datasets attached to the package. See examples on how to use Load() for details.\cr\cr

The data is provided through a variable named 'sampleMap' and is structured as expected from the 'Load()' function in the 's2dverification' package if was called as follows:\cr\cr
  \preformatted{
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
                    file_path = file.path('$STORE_FREQ$_mean', 
                      '$VAR_NAME$_3hourly/$VAR_NAME$_$START_DATE$.nc')
c <- ConfigAddEntry(c, 'observations', dataset_name = 'observation', 
                    var_name = 'tos', main_path = obs_data_path,
                    file_path = file.path('$STORE_FREQ$_mean',
                      '$VAR_NAME$/$VAR_NAME$_$YEAR$$MONTH$.nc')
ConfigFileSave(c, configfile, confirm = FALSE)

# Now we are ready to use Load().
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- Load('tos', c('experiment'), c('observation'), startDates, 
                   output = 'lonlat', latmin = 27, latmax = 48, lonmin = -12,
                   lonmax = 40, configfile = configfile)
  }
Check the documentation on 'Load()' in the package 's2dverification' for more information.
}
\usage{
data(sampleMap)
}
\format{
The data set provides with a variable named 'sampleMap'.\cr\cr

sampleMap$mod is an array that contains the experimental data and the dimension meanings and values are:\cr
  c(# of experimental datasets, # of members, # of starting dates, # of lead-times, # of latitudes, # of longitudes)\cr
  c(1, 3, 5, 60, 2, 3)\cr\cr

sampleMap$obs is an array that contains the observational data and the dimension meanings and values are:\cr
    c(# of observational datasets, # of members, # of starting dates, # of lead-times, # of latitudes, # of longitudes)\cr
    c(1, 1, 5, 60, 2, 3)\cr\cr

  sampleMap$lat is an array with the 2 latitudes covered by the data (see examples on Load() for details on why such low resolution).\cr\cr

  sampleMap$lon is an array with the 3 longitudes covered by the data (see examples on Load() for details on why such low resolution).
}
\keyword{datasets}

