#'Computes Probabilistic Information of a Forecast Relative to a Threshold or a Quantile
#'
#'Compute probabilistic bins of a set of forecast years ('fcyr') relative to 
#'the forecast climatology over the whole period of anomalies, optionally excluding 
#'the selected forecast years ('fcyr') or the forecast year for which the 
#'probabilistic bins are being computed (see 'compPeriod').
#'
#'@param ano Array of anomalies from Ano().\cr
#'  Must be of dimension (nexp/nobs, nmemb, nsdates, nleadtime, nlat,  nlon)
#'@param fcyr Indices of the forecast years of the anomalies which to compute 
#'  the probabilistic bins for, or 'all' to compute the bins for all the 
#'  years.\cr
#'  E.g., c(1:5), c(1, 4), 4 or 'all'.
#'@param thr Values used as thresholds to bin the anomalies.
#'@param quantile If quantile is TRUE (default), the threshold ('thr') 
#'  are quantiles.\cr
#'  If quantile is FALSE the thresholds ('thr') introduced are the absolute 
#'  thresholds of the bins.
#'@param posdates Position of the dimension in \code{ano} that corresponds to 
#'  the start dates (default = 3).
#'@param posdim Position of the dimension in \code{ano} which will be combined 
#'  with 'posdates' to compute the quantiles (default = 2, ensemble members).
#'@param compPeriod Three options: 
#'  "Full period"/"Without fcyr"/"Cross-validation" (The probabilities are 
#'  computed with the terciles based on ano/ano with all 'fcyr's 
#'  removed/cross-validation). The default is "Full period".
#'
#'@return Array with probabilistic information and dimensions:\cr
#'  c(length('thr') + 1, length(fcyr), nmemb/nparam, nmod/nexp/nobs, 
#'  nltime, nlat, nlon)\cr
#'  The values along the first dimension take values 0 or 1 depending on which 
#'  of the 'thr'+1 cathegories the forecast/observation at the corresponding 
#'  grid point, time step, member and starting date belongs to.
#'
#'@keywords datagen
#'@author History:\cr
#'1.0  -  2013  (F.Lienert)  -  Original code\cr
#'2.0  -  2014-03  (N. Gonzalez and V. Torralba, \email{veronica.torralba@@bsc.es})  -  Debugging
#'2.1  -  2017-02  (V. Torralba and N. Manubens, \email{veronica.torralba@@bsc.es})  -  Fix bug with cross-validation
#'@examples
#'# See examples on Load() to understand the first lines in this example
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
#'                   output = 'lonlat', latmin = 27, latmax = 48, 
#'                   lonmin = -12, lonmax = 40)
#'  }
#'  \dontshow{
#'startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
#'sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'),
#'                                                c('observation'), startDates,
#'                                                output = 'lonlat',
#'                                                latmin = 27, latmax = 48,
#'                                                lonmin = -12, lonmax = 40)
#'  }
#'clim <- Clim(sampleMap$mod, sampleMap$obs)
#'ano_exp <- Ano(sampleMap$mod, clim$clim_exp)
#'PB <- ProbBins(ano_exp, fcyr = 3, thr = c(1/3, 2/3), quantile = TRUE, posdates = 3,
#'               posdim = 2)
#'
#'@export
ProbBins <- function(ano, fcyr = 'all', thr, quantile = TRUE, posdates = 3,
                     posdim = 2, compPeriod = "Full period") {
  # define dimensions
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  nbdim <- length(dim(ano))

  if (nbdim < 7){
    ano <- Enlarge(ano, 7)
  }

  #remove posdates of posdim and supress duplicates
  posdim <- setdiff(posdim, posdates)
  nbpos <- length(posdim)
  #permute dimensions in ano
  if (posdates != 1 || posdim != 2) {
    dimnames_backup <- names(dim(ano))
    perm <- c(posdates, posdim, (1:7)[-c(posdates, posdim)])
    ano <- aperm(ano, perm)
    names(dim(ano)) <- dimnames_backup[perm]
    posdates <- 1
    posdim <- 2
  }
  dimano <- dim(ano)
  nsdates <- dimano[1]
  #calculate the number of elements on posdim dimension in ano
  nmemb <- 1
  if (nbpos > 0){
    for (idim in 2:(nbpos+1)){
      nmemb <- nmemb*dimano[idim]
    }
  }

  if (length(fcyr) == 1) {
    if (fcyr == 'all') {
      fcyr <- 1:nsdates
    }
  }
  nfcyr <- length(fcyr)

  if (compPeriod == "Cross-validation") {
    result <- NULL
    for (iyr in fcyr) {
      if (is.null(result)) {
        result <- ProbBins(ano, iyr, thr, quantile, 
                           posdates, posdim, "Without fcyr")
      } else {
        dimnames_backup <- names(dim(result))
        result <- abind(result, ProbBins(ano, iyr, thr, quantile, 
                                         posdates, posdim, "Without fcyr"), 
                        along = 2)
        names(dim(result)) <- dimnames_backup
      }
    }
    return(result)
  } else if (compPeriod %in% c("Full period", "Without fcyr")) {
    # separate forecast and hindcast
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    fore <- array(ano[fcyr, , , , , , ], dim = c(nfcyr,
                                                 dimano[2:7]))
    # the members and startdates are combined in one dimension
    sample_fore <- array(fore, dim=c(nfcyr*nmemb, dimano[(nbpos+2):7]))
    
    if(compPeriod == "Full period") {
      hind <- ano
      sample <- array(hind, dim=c(nsdates*nmemb, dimano[(nbpos+2):7]))
    } else if (compPeriod == "Without fcyr") {
      hind <- array(ano[-fcyr, , , , , , ], dim = c(nsdates-nfcyr,
                                                   dimano[2:7]))
      sample <- array(hind, dim=c((nsdates-nfcyr)*nmemb, dimano[(nbpos+2):7]))
    }
  
    #quantiles for each grid point and experiment
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if (quantile==TRUE){
      qum <- apply(sample, seq(2,7-nbpos,1), FUN=quantile,probs=thr,na.rm=TRUE,names=FALSE,type=8)                                        
    }else{
      qum<-array(thr,dim=c(length(thr), dimano[(nbpos+2):7]))
    }
  
    # This function assign the values to a category which is limited by the thresholds
    # It provides binary information
    
    counts <- function (dat, nbthr){
      thr <- dat[1:nbthr]
      data <-  dat[nbthr+1:(length(dat)-nbthr)]
      prob <- array(NA, dim=c(nbthr+1, length(dat)-nbthr))
      prob[1,]=1*(data <= thr[1])
      if(nbthr!=1){
        for (ithr in 2:(nbthr)){
          prob[ithr,]=1*((data > thr[ithr - 1]) & (data <= thr[ithr]))
        }
      }
      prob[nbthr+1,]=1*(data > thr[nbthr])
      return(prob) 
    }
        
    # The thresholds and anomalies are combined to use apply
    data <- abind(qum, sample_fore, along = 1)
    
    # PBF:Probabilistic bins of a forecast.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This array contains zeros and ones that indicate the category where your forecast is. 
  
    PBF <- array(apply(data, seq(2,7-nbpos,1), FUN=counts, nbthr=length(thr)),
                 dim=c(length(thr)+1, nfcyr, nmemb, dimano[(nbpos+2):nbdim]))
    
    names(dim(PBF)) <- c('bin', 'sdate', 'member', names(dim(ano))[(nbpos+2):nbdim])
    return(PBF)
  } else {
    stop("Parameter 'compPeriod' must be one of 'Full period', 'Without fcyr' or 'Cross-validation'.")
  }
}
