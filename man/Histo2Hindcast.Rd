\name{Histo2Hindcast}
\alias{Histo2Hindcast}
\title{
Chunks Long Simulations For Comparison With Hindcasts
}
\description{
This function reorganizes a long run (historical typically) with only one start date into chunks corresponding to a set of start dates. The expected input structure is the one output from \code{Load()} with 4 to 7 dimensions.
}
\usage{
Histo2Hindcast(varin, sdatesin, sdatesout, nleadtimesout)
}
\arguments{
  \item{varin}{
Input model or observational data:\cr
  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltimes) up to\cr
  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltimes, nlevel, nlat, nlon)
  }
  \item{sdatesin}{
Start date of the input matrix 'YYYYMMDD'.
  }
  \item{sdatesout}{
List of start dates of the output matrix c('YYYYMMDD', 'YYYYMMDD', ...).
  }
  \item{nleadtimesout}{
Number of leadtimes in the output matrix.
  }
}
\value{
A matrix with the same number of dimensions as the input one, the same dimensions 1 and 2 and potentially the same dimensions 5 to 7. Dimensions 3 and 4 are set by the arguments sdatesout and nleadtimesout.
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
startDates <- c('19901101')
sampleData <- Load('tos', c('experiment'), c('observation'), startDates,
                   leadtimemin = 1, leadtimemax = 4, output = 'areave',
                   latmin = 27, latmax = 48, lonmin = -12, lonmax = 40,
                   configfile = configfile)
  }
  \dontshow{
startDates <- c('19901101')
sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'),
                                                c('observation'), startDates,
                                                leadtimemin = 1,
                                                leadtimemax = 4,
                                                output = 'areave',
                                                latmin = 27, latmax = 48,
                                                lonmin = -12, lonmax = 40)
  }


start_dates_out <- c('19901101', '19911101', '19921101', '19931101', '19941101')
leadtimes_per_startdate <- 12
experimental_data <- Histo2Hindcast(sampleData$mod, startDates[1], 
                                    start_dates_out, leadtimes_per_startdate)
observational_data <- Histo2Hindcast(sampleData$obs, startDates[1], 
                                     start_dates_out, leadtimes_per_startdate)
PlotAno(experimental_data, observational_data, start_dates_out, 
        toptitle = paste('anomalies reorganized into shorter chunks'), 
        ytitle = 'K', fileout='tos_histo2hindcast.eps')
}
\author{
History:\cr
0.1  -  2012-11  (V. Guemas, \email{vguemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
