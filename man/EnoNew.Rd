\name{EnoNew}
\alias{EnoNew}
\title{Computes Effective Sample Size Following Guemas et al, BAMS, 2013b}
\description{
This function computes the effective number of independent values in the xdata array following the method described in Guemas V., Auger L., Doblas-Reyes F., JAMC, 2013. \code{EnoNew} provides similar functionality to \code{Eno} but with the added options to remove the linear trend or filter the frequency.}
\usage{
EnoNew(xdata, detrend = FALSE, filter = FALSE)
}
\arguments{
  \item{xdata}{a numeric vector}
  \item{detrend}{should the linear trend be removed from the data prior to the estimation of the equivalent number of independent values.}
  \item{filter}{should a filtering of the frequency peaks be applied prior to the estimation of the equivalent number of independant data.}
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
                   leadtimemin = 1, leadtimemax = 4, output = 'lonlat',
                   latmin = 27, latmax = 48, lonmin = -12, lonmax = 40,
                   configfile = configfile)
  }
  \dontshow{
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'),
                                                c('observation'), startDates,
                                                leadtimemin = 1,
                                                leadtimemax = 4,
                                                output = 'lonlat',
                                                latmin = 27, latmax = 48,
                                                lonmin = -12, lonmax = 40)
  }
eno <- EnoNew(sampleData$mod[1, 1, , 1, 2, 3])
print(eno)
}
\references{
Guemas V, Auger L, Doblas-Reyes FJ, Rust H, Ribes A, 2014, Dependencies in Statistical Hypothesis Tests for Climate Time Series. Bulletin of the American Meteorological Society, 95 (11), 1666-1667.
}
\author{
History:\cr
0.1  -  2012-06  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
