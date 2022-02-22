#'Computes Seasonal Means
#'
#'Computes seasonal means on timeseries organized in a array of any number of 
#'dimensions up to 10 dimensions where the time dimension is one of those 10 
#'dimensions.
#'
#'@param var Array containing the timeseries along one of its dimensions. 
#'@param posdim Dimension along which to compute seasonal means = Time 
#'  dimension.
#'@param monini an integer indicating the first month of the time series: 1 to 
#'  12.
#'@param moninf an integer indicating the month when to start the seasonal 
#'  means: 1 to 12.
#'@param monsup an integer indicating the month when to stop the seasonal 
#'  means: 1 to 12.
#'
#'@return Array with the same dimensions as var except along the posdim 
#'  dimension whose length corresponds to the number of seasons. Partial 
#'  seasons are not accounted for.
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-03  (V. Guemas, \email{virginie.guemas@bsc.es})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@bsc.es})  -  Formatting to CRAN
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
#'leadtimes_dimension <- 4
#'initial_month <- 11
#'mean_start_month <- 12
#'mean_stop_month <- 2
#'season_means_mod <- Season(sampleData$mod, leadtimes_dimension, initial_month,
#'                           mean_start_month, mean_stop_month)
#'season_means_obs <- Season(sampleData$obs, leadtimes_dimension, initial_month,
#'                           mean_start_month, mean_stop_month)
#'  \donttest{
#'PlotAno(season_means_mod, season_means_obs, startDates, 
#'        toptitle = paste('winter (DJF) temperatures'), ytitle = c('K'), 
#'        legends = 'ERSST', biglab = FALSE, fileout = 'tos_season_means.eps')
#'  }
#'@export
Season <- function(var, posdim = 4, monini, moninf, monsup) {
  while (monsup < moninf) {
    monsup <- monsup + 12
  }
  #
  #  Enlarge the size of var to 10 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
  }
  ntime <- dimsvar[posdim]
  enlvar <- Enlarge(var, 10)
  outdim <- c(dimsvar, array(1, dim = (10 - length(dimsvar))))
  u <- IniListDims(outdim, 10)
  v <- IniListDims(outdim, 10)
  #
  #  Compute the seasonal means 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  ind <- 1:ntime
  months <- ((ind - 1) + monini - 1) %% 12 + 1
  years <- ((ind - 1) + monini - 1) %/% 12

  for (jmon in moninf:monsup) {
    u[[posdim]] <- ind[which(months == ((jmon - 1) %% 12 + 1))]
    ind0 <- u[[posdim]][1]
    indf <- u[[posdim]][length(u[[posdim]])]
    if (indf > (ntime - (monsup - jmon))) {
      u[[posdim]] <- u[[posdim]][-which(u[[posdim]] == indf)]
    }
    if (ind0 < (jmon - moninf + 1)) {
      u[[posdim]] <- u[[posdim]][-which(u[[posdim]] == ind0)]
    } 
    if (jmon == moninf) { 
      nseas <- length(u[[posdim]])
      dimsvar[posdim] <- nseas
      outdim[posdim] <- nseas
      enlvarout <- array(0, dim = outdim)
    }
    v[[posdim]] <- 1:nseas
    enlvarout[v[[1]], v[[2]], v[[3]], v[[4]], v[[5]], v[[6]], v[[7]], v[[8]], 
              v[[9]], v[[10]]] <- enlvarout[v[[1]], v[[2]], v[[3]], v[[4]],
                                            v[[5]], v[[6]], v[[7]], v[[8]], 
                                            v[[9]], v[[10]]] + enlvar[u[[1]], 
                                  u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], 
                                  u[[7]], u[[8]], u[[9]], u[[10]]]
  }
  varout <- array(dim = dimsvar)
  varout[] <- enlvarout
  varout <- varout / (monsup - moninf + 1)
  #
  #  Outputs
  # ~~~~~~~~~
  #
  varout
}
