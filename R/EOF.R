#'Area-Weighted Empirical Orthogonal Function Analysis Using SVD
#'
#'Performs an area-weighted EOF analysis using SVD based on a covariance matrix 
#'by default, based on the correlation matrix if \code{corr} argument is set to 
#'\code{TRUE}.
#'
#'@param ano Array of anomalies with dimensions (number of timesteps, 
#'  number of latitudes, number of longitudes). NAs could exist but it should 
#'  be consistent along time_dim. That is, if one grid point has NAs, all the 
#'  time steps at this point should be NAs. 
#'@param lon Vector of longitudes of \code{ano}.
#'@param lat Vector of latitudes of \code{ano}.
#'@param neofs Number of modes to be kept. Default = 15.
#'@param corr Whether to base on a correlation matrix (\code{TRUE}) or on a 
#'  covariance matrix (default, \code{FALSE}).
#'
#'@return 
#'\item{EOFs}{
#'  An array of EOF patterns normalized to 1 (unitless) with dimensions 
#'  (number of modes, number of latitudes, number of longitues). 
#'  Multiplying \code{EOFs} by \code{PCs} gives the original reconstructed field.
#'}
#'\item{PCs}{
#'  An array of pincipal components with the units of the original field to 
#'  the power of 2, with dimensions (number of time steps, number of modes). 
#'  \code{PCs} contains already the percentage of explained variance so, 
#'  to reconstruct the original field it's only needed to multiply \code{EOFs} 
#'  by \code{PCs}.
#'} 
#'\item{var}{
#'  Percentage (%) of variance fraction of total variance explained by each 
#'  mode (number of modes).
#'}
#'\item{mask}{
#'  Mask with dimensions (number of latitudes, number of longitudes).
#'}
#'\item{wght}{
#'  Weights with dimensions (number of latitudes, number of longitudes).
#'}
#'
#'@seealso ProjectField, NAO, PlotBoxWhisker
#'@keywords datagen
#'@author History:\cr
#' 0.1  -  2012-10  (F. Lienert)  -  Original
#' code, inspired by R. Benestad's EOF() in R package clim.pact.\cr
#' 0.2  -  2014-03  (Lauriane Batte)  -  Bug-fixes:\cr
#'        1- Reversion of latitudes in the weights\cr
#'        2- Correlation matrix was used instead of covariance\cr
#'        3- Double use of the weights\cr
#'0.3  -  2014-03  (Virginie Guemas)  -  Bug-fixes:\cr
#'        1- Weight computation - division by sum of cos(lat)\cr
#'        2- Shuffling of EOFs in EOF.2 intermediate vector\cr
#'        3- Crash when neofs = 1 sorted out\cr
#'        4- Crash when neofs > nt sorted out\cr
#'0.4  -  2014-03  (Lauriane Batte)  -  Fixes:\cr
#'        1- BIG cleanup of code and clarification\cr
#'        2- Reduction of the number of transpositions and associated bug-fixes\cr
#'        4- Remove of the obsolete LINPACK options\cr
#'0.5  -  2014-04  (Virginie Guemas)  - Fixes:\cr
#'        1- Bug-fix in dimensions handling EOF composition restitutes now the
#'original field in all cases\cr
#'        2- Simplification of the convention transpose\cr
#'        3- Options to use the correlation matrix rather than the 
#'covariance matrix\cr
#'        4- Security checks\cr
#'        5- New normalization of PCs so that PC*EOF only reconstruct the 
#'original file\cr
#'        6- Weights = sqrt(cos(lat)) for ano so that covariance matrice 
#'weighted by cos(lat)\cr
#'        7- Division of EOF by weights so that the reconstruction is simply 
#'EOF * PC\cr
#'1.0  -  2016-03  (N. Manubens)  -  Formatting to R CRAN
#' 
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
#'# This example computes the EOFs along forecast horizons and plots the one that 
#'# explains the greatest amount of variability. The example data is very low 
#'# resolution so it does not make a lot of sense.
#'ano <- Ano_CrossValid(sampleData$mod, sampleData$obs)
#'eof <- EOF(Mean1Dim(ano$ano_exp, 2)[1, , 1, , ], sampleData$lon, sampleData$lat)
#'PlotEquiMap(eof$EOFs[1, , ], sampleData$lon, sampleData$lat)
#'
#'@importFrom stats sd
#'@export
EOF <- function(ano, lon, lat, neofs = 15, corr = FALSE) {
  # Checking ano
  if (!is.numeric(ano) || !is.array(ano)) {
    stop("Parameter 'ano' must be a numeric array.")
  }
  if (length(dim(ano)) != 3) {
    stop("'ano' must have dimensions c(n. start dates/forecast time steps/time samples, n. latitudes, n. longitudes).")
  }

  # Checking lon and lat
  if (!is.numeric(lon) || !is.numeric(lat)) {
    stop("'lon' and 'lat' must be numeric vectors.")
  }
  if (any(lon > 360 | lon < -360)) {
    .warning("some 'lon's out of the range [-360, 360].")
  }
  if (any(lat > 90 | lat < -90)) {
    stop("'lat' must contain values within the range [-90, 90].")
  }

  # Checking neofs
  if (!is.numeric(neofs)) {
    stop("'neofs' must be numeric.")
  }
  neofs <- round(neofs)

  # Checking corr
  if (!is.logical(corr)) {
    stop("'corr' must be either TRUE or FALSE.")
  }

  # Dimensions
  nlon <- length(lon)
  nlat <- length(lat)
  dim.dat <- dim(ano)
  ny <- dim.dat[2]
  nx <- dim.dat[3]
  nt <- dim.dat[1]
  
  # Security check
  if (ny != nlat) {
    stop("Inconsistent number of latitudes and input field dimensions.")
  }
  if (nx != nlon) {
    stop("Inconsistent number of longitudes and input field dimensions.")
  }
  
  # Check if all the time steps at one grid point are NA-consistent.
  # The grid point should have all NAs or no NA along time dim.
  if (any(is.na(ano))) {
    ano_latlon <- array(ano, dim = c(nt, ny * nx))  # [time, lat*lon]
    na_ind <- which(is.na(ano_latlon), arr.ind = T)
    if (dim(na_ind)[1] != nt * length(unique(na_ind[, 2]))) {
      stop("Detect certain grid points have NAs but not consistent across time ",
           "dimension. If the grid point is NA, it should have NA at all time step.")
    }
  }

  # Buildup of the mask
  mask <- ano[1, , ]
  mask[!is.finite(mask)] <- NA
  mask[is.finite(mask)] <- 1
  # Replace mask of NAs with 0s for EOF analysis.
  ano[!is.finite(ano)] <- 0
  
  # Area weighting. Weights for EOF; needed to compute the
  # fraction of variance explained by each EOFs
  wght <- array(cos(lat * pi/180), dim = (c(nlat, nlon)))
  
  # We want the covariance matrix to be weigthed by the grid
  # cell area so the anomaly field is weighted by its square
  # root since the covariance matrix equals transpose(ano)
  # times ano.
  wght <- sqrt(wght)
  ano <- ano * InsertDim(wght, 1, nt)
  
  # The use of the correlation matrix is done under the option
  # corr.
  if (corr == TRUE) {
    stdv <- apply(ano, sd, c(2, 3), na.rm = T)
    ano <- ano/InsertDim(stdv, 1, nt)
  }
  
  # Time/space matrix for SVD
  dim(ano) <- c(nt, ny * nx)
  dim.dat <- dim(ano)
  
  # 'transpose' means the array needs to be transposed before
  # calling La.svd for computational efficiency because the
  # spatial dimension is larger than the time dimension. This
  # goes with transposing the outputs of LA.svd also.
  if (dim.dat[2] > dim.dat[1]) {
    transpose <- TRUE
  } else {
    transpose <- FALSE
  }
  if (transpose) {
    pca <- La.svd(t(ano))
  } else {
    pca <- La.svd(ano)
  }
  
  # neofs is bounded
  neofs <- min(dim(ano), neofs)
  
  # La.svd conventions: decomposition X = U D t(V) La.svd$u
  # returns U La.svd$d returns diagonal values of D La.svd$v
  # returns t(V) !!  The usual convention is PC=U and EOF=V.
  # If La.svd is called for ano (transpose=FALSE case): EOFs:
  # $v PCs: $u If La.svd is called for t(ano) (transposed=TRUE
  # case): EOFs: t($u) PCs: t($v)
  
  if (transpose) {
    pca.EOFs <- t(pca$u)
    pca.PCs <- t(pca$v)
  } else {
    pca.EOFs <- pca$v
    pca.PCs <- pca$u
  }
  
  # The numbers of transposition is limited to neofs
  PC <- pca.PCs[, 1:neofs]
  EOF <- pca.EOFs[1:neofs, ]
  dim(EOF) <- c(neofs, ny, nx)
  
  # To sort out crash when neofs=1.
  if (neofs == 1) {
    PC <- InsertDim(PC, 2, 1)
  }
  
  # Computation of the % of variance associated with each mode
  W <- pca$d[1:neofs]
  tot.var <- sum(pca$d^2)
  var.eof <- 100 * pca$d[1:neofs]^2/tot.var
  
  for (e in 1:neofs) {
    
    # Factor to normalize the EOF.
    eof.patt.nn <- EOF[e, , ] * mask
    eof.patt.ms <- sum(eof.patt.nn^2, na.rm = TRUE)
    
    # Normalize the EOF
    eof.patt <- eof.patt.nn/eof.patt.ms
    
    # PC is multiplied by the normalization factor and the
    # weights, then the reconstruction is only EOF * PC (we have
    # multiplied ano by weight)
    eof.pc <- PC[, e] * eof.patt.ms * W[e]
    
    eof.patt <- eof.patt/wght
    
    EOF[e, , ] <- eof.patt
    PC[, e] <- eof.pc
  }
  return(list(EOFs = EOF, PCs = PC, var = var.eof, mask = mask, 
    wght = wght))
}
