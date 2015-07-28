\name{ACC}
\alias{ACC}
\title{
Computes Anomaly Correlation Coefficient (Spatial Correlation)
}
\description{
Matrix var_exp & var_obs should have dimensions (nexp/nobs, nsdates, nltimes, nlat, nlon) or (nexp/nobs, nsdates, nmember, nltimes, nlat, nlon).\cr
ACC computes the Anomaly Correlation Coefficient for the ensemble mean of each jexp in 1:nexp and each jobs in 1:nobs which gives nexp x nobs ACC for each startdate and each leadtime.\cr
A domain can be selected by providing the list of longitudes/latitudes (lon/lat) of the grid together with the corner of the domain:\cr
  lonlatbox = c(lonmin, lonmax, latmin, latmax)
}
\usage{
ACC(var_exp, var_obs, lon = NULL, lat = NULL, lonlatbox = NULL, conf = TRUE, 
    conftype = "parametric")
}
\arguments{
  \item{var_exp}{
Matrix of experimental anomalies with dimensions:\cr
  c(nexp, nsdates, nltimes, nlat, nlon)\cr
or\cr
  c(nexp, nsdates, nmembers, nltimes, nlat, nlon)\cr
  }
  \item{var_obs}{
Matrix of observational anomalies, same dimensions as var_exp except along the first dimension and the second if it corresponds to the member dimension.
  }
  \item{lon}{
Array of longitudes of the var_exp/var_obs grids, optional.
  }
  \item{lat}{
Array of latitudes of the var_exp/var_obs grids, optional.
  }
  \item{lonlatbox}{
Domain to select : c(lonmin, lonmax, latmin, latmax), optional.
  }
  \item{conf}{
TRUE/FALSE: confidence intervals or significance level provided or not.
  }
  \item{conftype}{
"parametric" provides a confidence interval for the ACC computed by a Fisher transformation and a significance level for the ACC from a one-sided student-T distribution.\cr
"bootstrap" provides a confidence interval for the ACC and MACC computed from bootstrapping on the members with 100 drawings with replacement.\cr
To guarantee the statistical robustness of the result, make sure that your experiments/oservations/startdates/leadtimes always have the same number of members.
  }
}
\value{
  \item{ACC}{
If conf set as TRUE, Matrix with dimensions:\cr
  c(nexp, nobs, nsdates, nleadtimes, 4) \cr
The fifth dimension of length 4 corresponds to the lower limit of the 95\% confidence interval, the ACC, the upper limit of the 95\% confidence interval and the 95\% significance level.\cr
If conf set as FALSE,  Anomaly Correlation Coefficient with dimensions:\cr
  c(nexp, nobs, nsdates, nleadtimes).  
  }
  \item{MACC}{
Mean Anomaly Correlation Coefficient with dimensions:\cr
  c(nexp, nobs, nleadtimes) 
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
sampleData$mod <- Season(sampleData$mod, 4, 11, 12, 2)
sampleData$obs <- Season(sampleData$obs, 4, 11, 12, 2)
clim <- Clim(sampleData$mod, sampleData$obs)
ano_exp <- Ano(sampleData$mod, clim$clim_exp)
ano_obs <- Ano(sampleData$obs, clim$clim_obs)
acc <- ACC(Mean1Dim(ano_exp, 2), Mean1Dim(ano_obs, 2))
PlotACC(acc$ACC, startDates)
}
\author{
History:\cr
0.1  -  2013-08  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN\cr
1.1  -  2013-09  (C. Prodhomme, \email{chloe.prodhomme at ic3.cat})  -  optimization\cr
1.2  -  2014-08  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Bug-fixes : handling of NA & selection of domain + Simplification of code\cr
1.3.0  -  2014-08  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Boostrapping over members\cr
1.3.1  -  2014-09  (C. Prodhomme, chloe.prodhomme at ic3.cat)  -  Add comments and minor style changes\cr
1.3.2  -  2015-02  (N. Manubens, nicolau.manubens at ic3.cat)  -  Fixed ACC documentation and examples\cr
}
\keyword{datagen}
