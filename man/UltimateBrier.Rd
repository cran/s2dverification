\name{UltimateBrier}
\alias{UltimateBrier}
\title{
Computes Brier Scores
}
\description{
Interface to compute probabilistic scores (Brier Score, Brier Skill Score) from data obtained from s2dverification.
}
\usage{
UltimateBrier(ano_exp, ano_obs, posdatasets = 1, posmemb = 2, posdates = 3, 
              quantile = TRUE, thr = c(5/100, 95/100), type = 'BS', 
              decomposition = TRUE)
}
\arguments{
  \item{ano_exp}{
Array of forecast anomalies, as provided by \code{Ano()}. Dimensions c(n. of experimental datasets, n. of members, n. of start dates, n. of forecast time steps, n. of latitudes, n. of longitudes). Dimensions in other orders are also supported. See parameters \code{posdatasets}, \code{posmemb} and \code{posdates}. 
  }
  \item{ano_obs}{
Array of observational reference anomalies, as provided by \code{Ano()}. Dimensions c(n. of observational reference datasets, n. of members, n. of start dates, n. of forecast time steps, n. of latitudes, n. of longitudes). Dimensions in other orders are also supported. See parameters \code{posdatasets}, \code{posmemb} and \code{posdates}.
  }
  \item{posdatasets}{
Expected position of dimension corresponding to the different evaluated datasets in input data (ano_exp and ano_obs). By default 1.
  }
  \item{posmemb}{
Expected position of dimension corresponding to members in input data (ano_exp and ano_obs). By default 2.
  }
  \item{posdates}{
Expected position of dimension corresponding to starting dates in input data (ano_exp and ano_obs). By default 3.
  }
  \item{quantile}{
Flag to stipulate whether a quantile (TRUE) or a threshold (FALSE) is used to estimate the forecast and observed probabilities. Takes TRUE by default.
  }
  \item{thr}{
Values to be used as quantiles if 'quantile' is TRUE or as thresholds if 'quantile' is FALSE. Takes by default \code{c(0.05, 0.95)} if 'quantile' is TRUE.
  }
  \item{type}{
Type of score desired. Can take the following values:
    \itemize{
      \item{'BS': Simple Brier Score.}
      \item{'FairEnsembleBS': Corrected Brier Score computed across ensemble members.}
      \item{'FairStartDatesBS': Corrected Brier Score computed across starting dates.}
      \item{'BSS': Simple Brier Skill Score.}
      \item{'FairEnsembleBSS': Corrected Brier Skill Score computed across ensemble members.}
      \item{'FairStartDatesBSS': Corrected Brier Skill Score computed across starting dates.}
    }
  }
  \item{decomposition}{
Flag to determine whether the decomposition of the Brier Score into its components should be provided (TRUE) or not (FALSE). Takes TRUE by default. The decomposition will be computed only if 'type' is 'BS' or 'FairStartDatesBS'.
  }
}
\value{
If 'type' is 'FairEnsembleBS', 'BSS', 'FairEnsemblesBSS' or 'FairStartDatesBSS' or 'decomposition' is FALSE and 'type' is 'BS' or 'FairStartDatesBS', the Brier Score or Brier Skill Score will be returned respectively.
If 'decomposition' is TRUE and 'type' is 'BS' or 'FairStartDatesBS' the returned value is a named list with the following entries:
  \itemize{
    \item{'BS': Brier Score.}
    \item{'REL': Reliability component.}
    \item{'UNC': Uncertainty component.}
    \item{'RES': Resolution component.}
  }
The dimensions of each of these arrays will be c(n. of experimental datasets, n. of observational reference datasets, n. of bins, the rest of input dimensions except for the ones pointed by 'posmemb' and 'posdates').
}
\examples{
# See ?Load for an explanation on the first part of this example.
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
                   leadtimemin = 1, leadtimemax = 4, output = 'lonlat',
                   latmin = 27, latmax = 48, lonmin = -12, lonmax = 40)
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
bs <- UltimateBrier(ano_exp, ano_obs)
bss <- UltimateBrier(ano_exp, ano_obs, type = 'BSS')
}
\author{
History:\cr
0.1  -  2015-05 (V. Guemas \email{virginie.guemas at bsc.es},\cr
                 C. Prodhomme \email{chloe.prodhomme at bsc.es},\cr
                 O. Bellprat \email{omar.bellprat at bsc.es}\cr
                 V. Torralba \email{veronica.torralba at bsc.es}\cr
                 N. Manubens, \email{nicolau.manubens at bsc.es})  -  First version
}
\keyword{datagen}
