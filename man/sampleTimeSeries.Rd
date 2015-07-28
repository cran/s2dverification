\name{sampleTimeSeries}
\docType{data}
\alias{sampleTimeSeries}
\title{Sample Of Observational And Experimental Data For Forecast Verification As Area Averages}
\description{
This data set provides area averaged data for the variable 'tos', i.e. sea surface temperature, over the mediterranean zone from the example datasets attached to the package. See examples on Load() for more details.\cr\cr

The data is provided through a variable named 'sampleTimeSeries' and is structured as expected from the 'Load()' function in the 's2dverification' package if was called as follows:\cr\cr
configfile <- paste0(tempdir(), '/sample.conf')\cr
ConfigFileCreate(configfile, confirm = FALSE)\cr
c <- ConfigFileOpen(configfile)\cr
c <- ConfigEditDefinition(c, 'DEFAULT_GRID', 'r20x10', confirm = FALSE)\cr
c <- ConfigEditDefinition(c, 'DEFAULT_VAR_MIN', '-1e19', confirm = FALSE)\cr
c <- ConfigEditDefinition(c, 'DEFAULT_VAR_MAX', '1e19', confirm = FALSE)\cr
c <- ConfigAddVar(c, '2d', 'tos')\cr
data_path <- system.file('sample_data', package = 's2dverification')\cr
exp_data_path <- paste0(data_path, '/model/$EXP_NAME$')\cr
obs_data_path <- paste0(data_path, '/$OBS_NAME$')\cr
c <- ConfigAddEntry(c, 'experiments', 'file-per-startdate',
     dataset_name = 'experiment', var_name = 'tos', main_path = exp_data_path,
     file_path = '$STORE_FREQ$_mean/$VAR_NAME$_3hourly/$VAR_NAME$_$START_DATE$.nc')
c <- ConfigAddEntry(c, 'observations', 'file-per-month',
     dataset_name = 'observation', var_name = 'tos', main_path = obs_data_path,
     file_path = '$STORE_FREQ$_mean/$VAR_NAME$/$VAR_NAME$_$YEAR$$MONTH$.nc')
ConfigFileSave(c, configfile, confirm = FALSE)\cr\cr

# Now we are ready to use Load().\cr
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')\cr
sampleData <- Load('tos', c('experiment'), c('observation'), startDates, \cr
                   output = 'areave', latmin = 27, latmax = 48, lonmin = -12,\cr
                   lonmax = 40, configfile = configfile)\cr\cr
Check the documentation on 'Load()' in the package 's2dverification' for more information.
}
\usage{
data(sampleTimeSeries)
}
\format{
The data set provides with a variable named 'sampleTimeSeries'.\cr\cr

sampleTimeSeries$mod is an array that contains the experimental data and the dimension meanings and values are:\cr
  c(# of experimental datasets, # of members, # of starting dates, # of lead-times)\cr
  c(1, 3, 5, 60)\cr\cr

sampleTimeSeries$obs is an array that contains the observational data and the dimension meanings and values are:\cr
  c(# of observational datasets, # of members, # of starting dates, # of lead-times)\cr
  c(1, 1, 5, 60)\cr\cr
    
sampleTimeSeries$lat is an array with the 2 latitudes covered by the data that was area averaged to calculate the time series (see examples on Load() for details on why such low resolution).\cr\cr

sampleTimeSeries$lon is an array with the 3 longitudes covered by the data that was area averaged to calculate the time series (see examples on Load() for details on why such low resolution).
}
\keyword{datasets}

