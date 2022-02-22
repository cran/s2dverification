#'Smoothes An Array Along A Dimension
#'
#'Smoothes an array of any number of dimensions along one of its dimensions.
#'
#'@param var Array to be smoothed along one of its dimension (typically the 
#'  forecast time dimension).
#'@param runmeanlen Running mean length in number of sampling units 
#'  (typically months).
#'@param numdimt Dimension to smooth.
#'
#'@return Array with same the dimensions as 'var' but smoothed along the 
#'  'numdimt'-th dimension.
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-03  (V. Guemas, \email{virginie.guemas@bsc.es})  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@bsc.es})  -  Formatting to R CRAN\cr
#'1.1  -  2015-05  (N. Manubens, \email{nicolau.manubens@bsc.es})  -  Adding
#'  security checks, fixing computation in cases where runmeanlen is odd and 
#'  making it able to work on arrays of any number of dimensions.
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
#'clim <- Clim(sampleData$mod, sampleData$obs)
#'ano_exp <- Ano(sampleData$mod, clim$clim_exp)
#'ano_obs <- Ano(sampleData$obs, clim$clim_obs)
#'runmean_months <- 12
#'dim_to_smooth <- 4  # Smooth along lead-times
#'smooth_ano_exp <- Smoothing(ano_exp, runmean_months, dim_to_smooth)
#'smooth_ano_obs <- Smoothing(ano_obs, runmean_months, dim_to_smooth)
#'  \donttest{
#'PlotAno(smooth_ano_exp, smooth_ano_obs, startDates, 
#'        toptitle = "Smoothed Mediterranean mean SST", ytitle = "K",
#'        fileout = "tos_smoothed_ano.eps")
#'  }
#'@import plyr
#'@export
Smoothing <- function(var, runmeanlen = 12, numdimt = 4) {
  # Check var
  if (!is.numeric(var)) {
    stop("Parameter 'var' must be a numeric vector or array.")
  }

  # Check runmeanlen
  if (!is.numeric(runmeanlen)) {
    stop("Parameter 'runmeanlen' must be numeric.")
  }
  runmeanlen <- round(runmeanlen)

  # Check numdimt
  if (!is.numeric(numdimt)) {
    stop("Parameter 'numdimt' must be numeric.")
  }
  numdimt <- round(numdimt)

  # Extra checks
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
    dim(var) <- dimsvar
  }
  if (numdimt > length(dimsvar) || numdimt < 1) {
    stop("Parameter 'numdimt' must be a number in the range [1, # of dimensions in 'var'].")
  }
  if (runmeanlen > dimsvar[numdimt] || runmeanlen < 1) {
    stop("Parameter 'runmeanlen' must be a number in the range [1, dim(var)[numdimt]].")
  }
  if (runmeanlen == dimsvar[numdimt] && runmeanlen %% 2 == 0) {
    stop("Parameter 'runmeanlen' is even and equal to the number of available data along dimension 'numdimt': no enough data. Either pick a narrower window or provide more data.")
  }

  # If only 1 dimension, add a dimension of lenght 1 so that code below 
  # still works
  flag_1dim <- FALSE
  if (length(dim(var)) == 1) {
    dim(var) <- c(1, dim(var))
    dimsvar <- dim(var)
    numdimt <- numdimt + 1
    flag_1dim <- TRUE
  }

  put <- function(x, along, indices, value) {
    nd <- length(dim(x))
    index <- as.list(rep(TRUE, nd))
    index[along] <- indices
    do.call("[<-", c(list(as.name("x")), index, list(value = as.name("value"))))
  }

  # If the window is of length > 1, we need to do computation. Else we just
  # return the original data
  if (runmeanlen == 1) {
    var
  } else {
    nmr1 <- floor(runmeanlen / 2)
    nltime <- dimsvar[numdimt]
    smoothed <- array(NA, dim = dimsvar)
    # We do a loop for all values along the time dimension for which we can 
    # take the comple window around them. Only for these values we have to do
    # computations. The other values are left to NA.
    # If the window size is even, we'll have to weight the half to the values
    # at the ends.
    if ((runmeanlen %% 2) == 0) {
      for (jtime in (1 + nmr1):(nltime - nmr1)) {
        # First, taking the two extreme values
        left <- take(var, numdimt, jtime - nmr1)
        right <- take(var, numdimt, jtime + nmr1)
        # Second, averaging the equally-weighted values around the centered
        mid <- take(var, numdimt, list((jtime - nmr1 + 1):(jtime + nmr1 - 1)))
        dims_mid <- dim(mid)
        dims_mid[numdimt] <- 1
        mid <- apply(mid, (1:length(dim(mid)))[-numdimt], function(x) sum(x) / runmeanlen)
        dim(mid) <- dims_mid
        # Finally, adding the center average to the ends weighted at half
        smoothed <- put(smoothed, numdimt, jtime, (left + right) / (runmeanlen * 2) + mid)
      }
    } else {
      for (jtime in (1 + nmr1):(nltime - nmr1)) {
        # Averaging the equally-weighted values around the centered
        mid <- take(var, numdimt, list((jtime - nmr1):(jtime + nmr1)))
        mid <- apply(mid, (1:length(dim(mid)))[-numdimt], function(x) sum(x) / runmeanlen)
        smoothed <- put(smoothed, numdimt, jtime, mid)
      }
    }
    if (flag_1dim) {
      dim(smoothed) <- dim(smoothed)[-1]
    }
    smoothed
  }
}
