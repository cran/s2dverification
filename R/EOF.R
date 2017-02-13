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
