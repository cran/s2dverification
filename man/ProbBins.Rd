% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProbBins.R
\name{ProbBins}
\alias{ProbBins}
\title{Computes Probabilistic Information of a Forecast Relative to a Threshold or a Quantile}
\usage{
ProbBins(
  ano,
  fcyr = "all",
  thr,
  quantile = TRUE,
  posdates = 3,
  posdim = 2,
  compPeriod = "Full period"
)
}
\arguments{
\item{ano}{Array of anomalies from Ano().\cr
Must be of dimension (nexp/nobs, nmemb, nsdates, nleadtime, nlat,  nlon)}

\item{fcyr}{Indices of the forecast years of the anomalies which to compute 
the probabilistic bins for, or 'all' to compute the bins for all the 
years.\cr
E.g., c(1:5), c(1, 4), 4 or 'all'.}

\item{thr}{Values used as thresholds to bin the anomalies.}

\item{quantile}{If quantile is TRUE (default), the threshold ('thr') 
are quantiles.\cr
If quantile is FALSE the thresholds ('thr') introduced are the absolute 
thresholds of the bins.}

\item{posdates}{Position of the dimension in \code{ano} that corresponds to 
the start dates (default = 3).}

\item{posdim}{Position of the dimension in \code{ano} which will be combined 
with 'posdates' to compute the quantiles (default = 2, ensemble members).}

\item{compPeriod}{Three options: 
"Full period"/"Without fcyr"/"Cross-validation" (The probabilities are 
computed with the terciles based on ano/ano with all 'fcyr's 
removed/cross-validation). The default is "Full period".}
}
\value{
Array with probabilistic information and dimensions:\cr
 c(length('thr') + 1, length(fcyr), nmemb/nparam, nmod/nexp/nobs, 
 nltime, nlat, nlon)\cr
 The values along the first dimension take values 0 or 1 depending on which 
 of the 'thr'+1 cathegories the forecast/observation at the corresponding 
 grid point, time step, member and starting date belongs to.
}
\description{
Compute probabilistic bins of a set of forecast years ('fcyr') relative to 
the forecast climatology over the whole period of anomalies, optionally excluding 
the selected forecast years ('fcyr') or the forecast year for which the 
probabilistic bins are being computed (see 'compPeriod').
}
\examples{
# See examples on Load() to understand the first lines in this example
 \dontrun{
data_path <- system.file('sample_data', package = 's2dverification')
expA <- list(name = 'experiment', path = file.path(data_path,
            'model/$EXP_NAME$/$STORE_FREQ$_mean/$VAR_NAME$_3hourly',
            '$VAR_NAME$_$START_DATE$.nc'))
obsX <- list(name = 'observation', path = file.path(data_path,
            '$OBS_NAME$/$STORE_FREQ$_mean/$VAR_NAME$',
            '$VAR_NAME$_$YEAR$$MONTH$.nc'))

# Now we are ready to use Load().
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- Load('tos', list(expA), list(obsX), startDates,
                  output = 'lonlat', latmin = 27, latmax = 48, 
                  lonmin = -12, lonmax = 40)
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
2.0  -  2014-03  (N. Gonzalez and V. Torralba)  -  Debugging
2.1  -  2017-02  (V. Torralba and N. Manubens)  -  Fix bug with cross-validation
}
\keyword{datagen}
