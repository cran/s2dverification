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
