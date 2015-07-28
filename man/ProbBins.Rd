\name{ProbBins}
\alias{ProbBins}
\title{
Computes probabilistic information of a forecast relative to a threshold or a quantile. 
}
\description{
Compute probabilistic bins of a forecast year (fcyr) relative to the forecast climatology over the whole period of anomalies excluding the selected forecast year (fcyr).
}
\usage{
ProbBins(ano, fcyr, thr, quantile=TRUE, posdates=3, posdim=2, 
         compPeriod = "Full period")
}
\arguments{
  \item{ano}{
Array of anomalies from Ano.\cr
Must be of dimension (nexp, nmemb, nsdates, nleadtime, nlon,  nlat)
  }
  \item{fcyr}{
Number of the forecast year of the anomalies.
  }
  \item{thr}{
values used to bin the anomalies.
  }
  \item{quantile}{
If quantile=TRUE (default), the threshold (thr) are quantiles.\cr
If quantile=FALSE the threshold (thr) introduced are the specified thresholds.
  }
  \item{posdates}{
Position of the associated dimension with the start dates (default=3).
  }
  \item{posdim}{
Position of the dimension which will be combined with posdates to compute the quantiles (default=2).
  }
  \item{compPeriod}{
Three options: "Full period"/"Without fcyr"/"cross-validation" (The probabilities are computed with the terciles based on ano/ano with the removed fcyr/cross-validation)(default=Full period).
  }
}
\value{
Matrix with probabilistic information and dimensions:\cr
  c(length(thr+1), nfcyr ,nmemb/nparam, nmod/nexp/nobs, nfcyr, nltime, nlat, nlon)
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
clim <- Clim(sampleMap$mod, sampleMap$obs)
ano_exp <- Ano(sampleMap$mod, clim$clim_exp)
PB <- ProbBins(ano_exp, fcyr = 3, thr = c(1/3, 2/3), quantile = TRUE, posdates = 3,
posdim = 2)
 
}
\author{
History:\cr
1.0  -  2013  (F.Lienert)  -  Original code\cr
2.0  -  2014-03  (N. Gonzalez and V.Torralba, \email{veronica.torralba at ic3.cat})  -  debugging
}
\keyword{datagen}
