#'Computes Brier Scores
#'
#'Interface to compute probabilistic scores (Brier Score, Brier Skill Score) 
#'from data obtained from s2dverification.
#'
#'@param ano_exp Array of forecast anomalies, as provided by \code{Ano()}. 
#'  Dimensions c(n. of experimental datasets, n. of members, n. of start dates, 
#'  n. of forecast time steps, n. of latitudes, n. of longitudes). Dimensions 
#'  in other orders are also supported. See parameters \code{posdatasets}, 
#'  \code{posmemb} and \code{posdates}.
#'@param ano_obs Array of observational reference anomalies, as provided by 
#'  \code{Ano()}. Dimensions c(n. of observational reference datasets, 
#'  n. of members, n. of start dates, n. of forecast time steps, 
#'  n. of latitudes, n. of longitudes). Dimensions in other orders are also 
#'  supported. See parameters \code{posdatasets}, \code{posmemb} and 
#'  \code{posdates}.
#'@param posdatasets Expected position of dimension corresponding to the 
#'  different evaluated datasets in input data (ano_exp and ano_obs). 
#'  By default 1.
#'@param posmemb Expected position of dimension corresponding to members in 
#'  input data (ano_exp and ano_obs). By default 2.
#'@param posdates Expected position of dimension corresponding to starting 
#'  dates in input data (ano_exp and ano_obs). By default 3.
#'@param quantile Flag to stipulate whether a quantile (TRUE) or a threshold 
#'  (FALSE) is used to estimate the forecast and observed probabilities. 
#'  Takes TRUE by default.
#'@param thr Values to be used as quantiles if 'quantile' is TRUE or as 
#'  thresholds if 'quantile' is FALSE. Takes by default \code{c(0.05, 0.95)} 
#'  if 'quantile' is TRUE.
#'@param type Type of score desired. Can take the following values:
#'\itemize{
#'  \item{'BS': Simple Brier Score.}
#'  \item{'FairEnsembleBS': Corrected Brier Score computed across ensemble 
#'    members.}
#'  \item{'FairStartDatesBS': Corrected Brier Score computed across starting 
#'    dates.}
#'  \item{'BSS': Simple Brier Skill Score.}
#'  \item{'FairEnsembleBSS': Corrected Brier Skill Score computed across 
#'    ensemble members.}
#'  \item{'FairStartDatesBSS': Corrected Brier Skill Score computed across 
#'    starting dates.}
#'}
#'@param decomposition Flag to determine whether the decomposition of the 
#'  Brier Score into its components should be provided (TRUE) or not (FALSE). 
#'  Takes TRUE by default. The decomposition will be computed only if 'type' 
#'  is 'BS' or 'FairStartDatesBS'.
#'@return 
#'If 'type' is 'FairEnsembleBS', 'BSS', 'FairEnsemblesBSS' or 
#''FairStartDatesBSS' or 'decomposition' is FALSE and 'type' is 'BS' or 
#''FairStartDatesBS', the Brier Score or Brier Skill Score will be returned 
#'respectively.
#'If 'decomposition' is TRUE and 'type' is 'BS' or 'FairStartDatesBS' the 
#'returned value is a named list with the following entries:
#'  \itemize{
#'    \item{'BS': Brier Score.}
#'    \item{'REL': Reliability component.}
#'    \item{'UNC': Uncertainty component.}
#'    \item{'RES': Resolution component.}
#'  }
#'The dimensions of each of these arrays will be c(n. of experimental datasets, 
#'n. of observational reference datasets, n. of bins, the rest of input 
#'dimensions except for the ones pointed by 'posmemb' and 'posdates').
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2015-05 (V. Guemas \email{virginie.guemas@@bsc.es},\cr
#'                 C. Prodhomme \email{chloe.prodhomme@@bsc.es},\cr
#'                 O. Bellprat \email{omar.bellprat@@bsc.es}\cr
#'                 V. Torralba \email{veronica.torralba@@bsc.es}\cr
#'                 N. Manubens, \email{nicolau.manubens@@bsc.es})  -  First version
#'@examples
#'# See ?Load for an explanation on the first part of this example.
#'  \dontrun{
#'data_path <- system.file('sample_data', package = 's2dverification')
#'expA <- list(name = 'experiment', path = file.path(data_path, 
#'             'model/$EXP_NAME$/$STORE_FREQ$_mean/$VAR_NAME$_3hourly',
#'             '$VAR_NAME$_$START_DATE$.nc'))
#'obsX <- list(name = 'observation', path = file.path(data_path,
#'             '$OBS_NAME$/$STORE_FREQ$_mean/$VAR_NAME$',
#'             '$VAR_NAME$_$YEAR$$MONTH$.nc'))
#'
#'# Now we are ready to use Load().
#'startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
#'sampleData <- Load('tos', list(expA), list(obsX), startDates, 
#'                   leadtimemin = 1, leadtimemax = 4, output = 'lonlat',
#'                   latmin = 27, latmax = 48, lonmin = -12, lonmax = 40)
#'  }
#'  \dontshow{
#'startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
#'sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'),
#'                                                c('observation'), startDates,
#'                                                leadtimemin = 1,
#'                                                leadtimemax = 4,
#'                                                output = 'lonlat',
#'                                                latmin = 27, latmax = 48,
#'                                                lonmin = -12, lonmax = 40)
#'  }
#'sampleData$mod <- Season(sampleData$mod, 4, 11, 12, 2)
#'sampleData$obs <- Season(sampleData$obs, 4, 11, 12, 2)
#'clim <- Clim(sampleData$mod, sampleData$obs)
#'ano_exp <- Ano(sampleData$mod, clim$clim_exp)
#'ano_obs <- Ano(sampleData$obs, clim$clim_obs)
#'bs <- UltimateBrier(ano_exp, ano_obs)
#'bss <- UltimateBrier(ano_exp, ano_obs, type = 'BSS')
#'@import SpecsVerification plyr
#'@export
UltimateBrier <- function(ano_exp, ano_obs, posdatasets = 1, posmemb = 2, 
                          posdates = 3, quantile = TRUE, 
                          thr = c(5/100, 95/100), type = 'BS', 
                          decomposition = TRUE) {
  # Checking ano_exp
  if (!is.numeric(ano_exp) || !is.array(ano_exp)) {
    stop("Parameter 'ano_exp' must be a numeric array.")
  }
  if (length(dim(ano_exp)) < 3) {
    stop("'ano_exp' must have at least the dimensions c(n. experimental data sets, n. members, n. start dates/forecast time steps/time steps).")
  }

  # Checking ano_obs
  if (!is.numeric(ano_obs) || !is.array(ano_obs)) {
    stop("Parameter 'ano_obs' must be a numeric array.")
  }
  if (length(dim(ano_obs)) < 3) {
    stop("'ano_obs' must have at least the dimensions c(n. observational data sets, n. obs. members, n. start dates/forecast time steps/time steps).")
  }

  # Checking consistency in ano_exp and ano_obs
  if (!(length(dim(ano_exp)) == length(dim(ano_obs)))) {
    stop("'ano_obs' and 'ano_exp' must have the same number of dimensions.")
  }
  if (!identical(dim(ano_exp)[c(1:length(dim(ano_exp)))[-c(posdatasets, posmemb)]], 
                 dim(ano_obs)[c(1:length(dim(ano_obs)))[-c(posdatasets, posmemb)]])) {
    stop("'ano_obs' and 'ano_exp' must have dimensions of equal size except for the datasets and members dimensions.")
  }
  if (dim(ano_obs)[posmemb] != 1) {
    stop("Only observations with one member are supported. dim(ano_obs)[posmemb] should be 1")
  }
  if (any(is.na(ano_exp)) || any(is.na(ano_obs))) {
    stop("There can't be NA values in 'ano_exp' and 'ano_obs'")
  }

  # Checking posdatasets
  if (!is.numeric(posdatasets)) {
    stop("Parameter 'posdatasets' must be an integer.")
  } else {
    posdatasets <- round(posdatasets)
    if (posdatasets > length(dim(ano_exp))) {
      stop("Parameter 'posdatasets' exceeds the number of dimensions of the provided anomalies.")
    }
  }

  # Checking posmemb
  if (!is.numeric(posmemb)) {
    stop("Parameter 'posmemb' must be an integer.")
  } else {
    posmemb <- round(posmemb)
    if (posmemb > length(dim(ano_exp))) {
      stop("Parameter 'posmemb' exceeds the number of dimensions of the provided anomalies.")
    }
  }

  # Checking posdates
  if (!is.numeric(posdates)) {
    stop("Parameter 'posdates' must be an integer.")
  } else {
    posdates <- round(posdates)
    if (posdates > length(dim(ano_exp))) {
      stop("Parameter 'posdates' exceeds the number of dimensions of the provided anomalies.")
    }
  }

  if (posdatasets == posmemb || posdatasets == posdates || posmemb == posdates) {
    stop("Parameters 'posdatasets', 'posmemb' and 'posdates' must all point to different dimensions.")
  }

  # Checking quantile
  if (!is.logical(quantile)) {
    stop("Parameter 'quantile' must be either TRUE or FALSE.")
  }

  # Checking thr
  if (!is.numeric(thr)) {
    stop("Parameter 'thr' must be a numerical vector.")
  }
  if (quantile) {
    if (!all(thr <= 1 & thr >= 0)) {
      stop("All quantiles specified in parameter 'thr' must fall in the range [0, 1].")
    }
  }

  # Checking type
  if (!(type %in% c('BS', 'FairEnsembleBS', 'FairStartDatesBS', 'BSS', 'FairEnsembleBSS', 'FairStartDatesBSS'))) {
    stop("'type' not supported")
  }

  # Checking decomposition
  if (!is.logical(decomposition)) {
    stop("Parameter 'decomposition' must be either TRUE or FALSE.")
  }

  #library(SpecsVerification)
  #library(plyr)

  ## The functions used to compute FEBSS and FEBS have as input anomalies from Ano()
  ## dims: c(datasets, members, startdates, leadtimes, latitudes, longitudes)
  if (type %in% c('FairEnsembleBSS', 'FairEnsembleBS')) {
    input_exp <- ano_exp
    input_obs <- ano_obs
    input_posdatasets <- posdatasets
    input_posmemb <- posmemb
    input_posdates <- posdates
  ## The functions used to compute the other scores receive data from ProbBins()
  ## dims: c(bins, startdates, members, datasets, leadtimes, latitudes, longitudes)
  } else {
    # Calculate probabilities of data with ProbBins and average over members
    exp_probs <- array(Mean1Dim(ProbBins(ano_exp, 1:dim(ano_exp)[posdates], thr, quantile, posdates, posmemb), 3),
                       dim = c(length(thr) + 1, dim(ano_exp)[posdates], 1, dim(ano_exp)[-c(posmemb, posdates)]))
    obs_probs <- array(Mean1Dim(ProbBins(ano_obs, 1:dim(ano_obs)[posdates], thr, quantile, posdates, posmemb), 3),
                       dim = c(length(thr) + 1, dim(ano_obs)[posdates], 1, dim(ano_obs)[-c(posmemb, posdates)]))
    input_exp <- exp_probs
    input_obs <- obs_probs
    input_posdatasets <- 4
    input_posmemb <- 3
    input_posdates <- 2
  }

  ## Here we define the function 'f' for each type of score (read further for more info).
  if (type == 'FairEnsembleBSS') {
    size_ens_ref <- prod(dim(ano_obs)[c(posmemb, posdates)])
    f <- function(x) {
      ens_ref <- matrix(do.call("[", c(list(x = x), indices_obs)), size_ens_ref, size_ens_ref, byrow = TRUE)
      sapply(c(thr, 1), function(tau) {
        FairBrierSs(t(do.call("[", c(list(x = x), indices_exp))) > tau, 
                    ens_ref > tau, 
                    do.call("[", c(list(x = x), indices_obs)) > tau)['skillscore']
      })
    }
  } else {
    if (type == 'FairEnsembleBS') {
      f <- function(x) sapply(c(thr, 1), function(tau) mean(FairBrier(t(do.call("[", c(list(x = x), indices_exp))) > tau, do.call("[", c(list(x = x), indices_obs)) > tau), na.rm = TRUE))
    } else {
      if (type == 'BS') {
        f <- function(x) as.vector(BrierScoreDecomposition(do.call("[", c(list(x = x), indices_exp)), do.call("[", c(list(x = x), indices_obs)))[1, ])
      } else if (type == 'FairStartDatesBS') {
        f <- function(x) unlist(BrierScore(do.call("[", c(list(x = x), indices_obs)), do.call("[", c(list(x = x), indices_exp)))[c('rel', 'res', 'unc')], use.names = FALSE)
      } else if (type == 'BSS') {
        f <- function(x) BrierScore(do.call("[", c(list(x = x), indices_obs)), do.call("[", c(list(x = x), indices_exp)))$bss_res
      } else if (type == 'FairStartDatesBSS') {
        f <- function(x) BrierScore(do.call("[", c(list(x = x), indices_obs)), do.call("[", c(list(x = x), indices_exp)))$bss_gres
      }
    }
  }

  ## We will calculate score for each exp, obs, bin, leadtime, latitude and longitude
  ## So we create array to store results
  ## If we calculate a BS we will store its rel, res and unc
  if (type %in% c('BS', 'FairStartDatesBS')) {
    result <- array(dim = c(dim(ano_exp)[posdatasets], dim(ano_obs)[posdatasets], 3, length(thr) + 1, dim(ano_exp)[-c(posdatasets, posmemb, posdates)]))
  } else {
    result <- array(dim = c(dim(ano_exp)[posdatasets], dim(ano_obs)[posdatasets], length(thr) + 1, dim(ano_exp)[-c(posdatasets, posmemb, posdates)]))
  }
  ## In order to be able to use apply, we put data of each exp and obs in a single array,
  ## all merged over the members dimension.
  indices_exp <- as.list(rep(TRUE, length(dim(input_exp))))
  indices_exp[[input_posmemb]] <- c(1:dim(input_exp)[input_posmemb])
  indices_exp <- indices_exp[c(input_posdatasets, input_posmemb, input_posdates)]
  indices_obs <- as.list(rep(TRUE, length(dim(input_obs))))
  indices_obs[[input_posmemb]] <- c(1:dim(input_obs)[input_posmemb]) + dim(input_exp)[input_posmemb]
  indices_obs <- indices_obs[c(input_posdatasets, input_posmemb, input_posdates)]
  out_indices <- as.list(rep(TRUE, length(dim(result))))
  for (n_obs in 1:dim(ano_obs)[posdatasets]) {
    for (n_exp in 1:dim(ano_exp)[posdatasets]) {
      data <- abind(take(input_exp, input_posdatasets, n_exp), 
                    take(input_obs, input_posdatasets, n_obs), 
                    along = input_posmemb)
      out_indices[c(1, 2)] <- c(n_exp, n_obs)
      ## We apply function 'f' to data of each couple of exp and obs, merged in a single array.
      ## This data will have dimensions
      ##   (1, nmembexp + nmembobs, nsdates, nltimes, nlat, nlon)
      ## or
      ##   (nbins, nsdates, nmembexp + nmembobs, 1, nltimes, nlat, nlon)
      ## depending on the input type.
      ## The function 'f' is applied along all dimensions but (datasets, members and sdates)
      ## so the produced output by apply is at least of dimensions
      ##   (nltimes, nlat, nlon)
      ## or
      ##   (nbins, nltimes, nlat, nlon)
      ## depending on the input type (FairEnsembleBS and FairEnsembleBSS 
      ## will have at lest the first set of dimensions, 
      ## other scores will have at least the second).
      ## So 'f' must have as input an array of dims (1, nmembexp + nmembobs, nsdates).
      ## 'indices_exp' and 'indices_obs' will pick for us the input data corresponding to
      ## experiments or observations respectively.
      ## In order to match with dimensions of 'result', the apply() must have as
      ## output an array of dims (nbins, nltimes, nlat, nlon)
      ##                      or (3, nbins, nltimes, nlat, nlon) if calculating BS or FairStartDatesBS.
      ## All in all, after looking at apply()'s 'at least' output
      ## dimensions and at apply()'s required output dimensions:
      ##   'f' must have as output a vector of length nbins for FEBS and FEBSS
      ##   'f' must have as output a vector of length 3 (corresponding to rel, res and unc) for BS and FSDBS
      ##   'f' must have as output a single numeric element for other scores
      result <- do.call("[<-", c(list(x = result), out_indices, list(value = apply(data, c(1:length(dim(data)))[-c(input_posdatasets, input_posmemb, input_posdates)], f))))
    }
  }

  if (type %in% c('BSS', 'FairStartDatesBSS', 'FairEnsembleBSS')) {
    result
  } else {
    if (decomposition && type != 'FairEnsembleBS') {
      rel <- take(result, 3, 1)
      dim(rel) <- dim(rel)[-3]
      res <- take(result, 3, 2)
      dim(res) <- dim(res)[-3]
      unc <- take(result, 3, 3)
      dim(unc) <- dim(unc)[-3]
      bs <- rel - res + unc
      list(BS = bs, REL = rel, UNC = unc, RES = res)
    } else {
      result <- take(result, 3, 1) - take(result, 3, 2) + take(result, 3, 3)
      dim(result) <- dim(result)[-3]
      result
    }
  }
}

