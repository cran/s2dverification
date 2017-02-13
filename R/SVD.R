SVD <- function(vary, varx, laty = NULL, latx = NULL, nmodes = 15,
                corr = FALSE, weight = TRUE) {
  # Checking vary
  if (!is.numeric(vary) || !is.array(vary)) {
    stop("Parameter 'vary' must be a numeric array.")
  }
  if (length(dim(vary)) != 3) {
    stop("'vary' must have dimensions c(n. start dates/forecast time steps/time samples, n. latitudes, n. longitudes).")
  }

  # Checking varx
  if (!is.numeric(varx) || !is.array(varx)) {
    stop("Parameter 'varx' must be a numeric array.")
  }
  if (length(dim(varx)) != 3) {
    stop("'varx' must have dimensions c(n. start dates/forecast time steps/time samples, n. latitudes, n. longitudes).")
  }

  # Checking weight
  if (!is.logical(weight)) {
    stop("'weight' must be either TRUE or FALSE.")
  }

  if (weight) {
    # Checking laty and latx
    if (!is.numeric(laty) || !is.numeric(latx)) {
      stop("'laty' and 'latx' must be numeric vectors.")
    }
    if (any(laty > 90 | laty < -90)) {
      stop("'laty' must contain values within the range [-90, 90].")
    }
    if (any(latx > 90 | latx < -90)) {
      stop("'latx' must contain values within the range [-90, 90].")
    }
  }

  # Checking nmodes
  if (!is.numeric(nmodes)) {
    stop("'nmodes' must be numeric.")
  }
  nmodes <- round(nmodes)

  # Checking corr
  if (!is.logical(corr)) {
    stop("'corr' must be either TRUE or FALSE.")
  }

  dy <- dim(vary)
  #if (length(dy) == 3) {
    nlony <- dy[3]
    nlaty <- dy[2]
    nsy <- nlony * nlaty
    nt <- dy[1]
  #} else {
  #  nlony <- 1
  #  nlaty <- dy[2]
  #  nsy <- nlony * nlaty
  #  nt <- dy[1]
  #}

  # Security check
  if (weight) {
    if (nlaty != length(laty)) {
      stop("Inconsistent number of latitudes ('laty') and input field dimensions ('vary').")
    }
    # Area weighting. Weights for EOF; needed to compute the
    # fraction of variance explained by each EOFs
    wghty <- array(cos(laty * pi/180), dim = (c(nlaty, nlony)))
    
    # We want the covariance matrix to be weigthed by the grid
    # cell area so the anomaly field is weighted by its square
    # root since the covariance matrix equals transpose(ano)
    # times ano.
    wghty <- sqrt(wghty)
    vary <- vary * InsertDim(wghty, 1, nt)
  }

  # The use of the correlation matrix is done under the option
  # corr.
  if (corr == TRUE) {
    stdv <- apply(vary, sd, c(2, 3), na.rm = TRUE)
    vary <- vary/InsertDim(stdv, 1, nt)
  }

  Y <- vary # time, lat, lon
  original_dims <- dim(Y)
  dim(Y) <- c(nt, nsy)
  y_na <- apply(Y, 2, function(time_series) if (any(is.na(time_series))) TRUE else FALSE)
  Yyes <- Y[, which(!y_na)] # time, space_reduced
  dim(Y) <- original_dims

  dx <- dim(varx)
  #if (length(dx) == 3) {
    nlonx <- dx[3]
    nlatx <- dx[2]
    nsx <- nlonx * nlatx
    nt2 <- dx[1]
  #} else {
  #  nlonx <- 1
  #  nlatx <- dx[2]
  #  nsx <- nlonx * nlatx
  #  nt2 <- dx[1]
  #}

  # Security check
  if (nt != nt2) {
    stop("Inconsistent number of time steps in input fields 'vary' and 'varx'.")
  }

  if (weight) {
    if (nlatx != length(latx)) {
      stop("Inconsistent number of latitudes ('latx') and input field dimensions ('varx').")
    }
    # Area weighting. Weights for EOF; needed to compute the
    # fraction of variance explained by each EOFs
    wghtx <- array(cos(latx * pi/180), dim = (c(nlatx, nlonx)))
    
    # We want the covariance matrix to be weigthed by the grid
    # cell area so the anomaly field is weighted by its square
    # root since the covariance matrix equals transpose(ano)
    # times ano.
    wghtx <- sqrt(wghtx)
    varx <- varx * InsertDim(wghtx, 1, nt)
  }

  # The use of the correlation matrix is done under the option
  # corr.
  if (corr == TRUE) {
    stdv <- apply(varx, sd, c(2, 3), na.rm = TRUE)
    varx <- varx/InsertDim(stdv, 1, nt)
  }

  X <- varx  # ntime, nlat, nlon
  original_dims <- dim(X)
  dim(X) <- c(nt, nsx)
  x_na <- apply(X, 2, function(time_series) if (any(is.na(time_series))) TRUE else FALSE)
  Xyes <- X[, which(!x_na)]  # ntime, nspace_reduced
  dim(X) <- original_dims

  # covariance matrix
  
  covXY <- (t(Xyes) %*% Yyes)/(nt - 1)
  
  # MCA analysis
  
  MCA <- svd(covXY)
  
  # outputs
 
  SV <- MCA$d  #singular values
  SV2 <- SV * SV
  S <- sum(SV2)
  SCF <- SV2/S  #squared covariance fraction
  nm <- min(nmodes, length(SV))  #number of modes to be retained
  SCF <- SCF[1:nm] * 100
  
  Ux <- MCA$u
  Vy <- MCA$v
  Un <- Ux[, 1:nm]
  Vn <- Vy[, 1:nm]
  ECu <- t(Xyes %*% Un)
  ECv <- t(Yyes %*% Vn)  #(nm, nt)
  ECu_std <- array(dim = c(nm, nt))
  ECv_std <- array(dim = c(nm, nt))
  for (i in 1:nm) {
    Su <- sd(ECu[i, ])
    ECu_std[i, ] <- ECu[i, ]/Su
    Sv <- sd(ECv[i, ])
    ECv_std[i, ] <- ECv[i, ]/Sv
  }
  
  RUV <- array(dim = nm)
  for (i in 1:nm) {
    RUV[i] <- cor(ECu_std[i, ], ECv_std[i, ])
  }
 
  MCAu <- array(dim = c(nm, nlatx, nlonx))
  MCAv <- array(dim = c(nm, nlaty, nlony))
  for (m in 1:nm) {
    MCAu[m, , ] <- apply(X, c(2, 3), 
      function(time_series) {
        if (any(is.na(time_series))) {
          NA
        } else {
          lm(time_series ~ ECu_std[m, ])$coefficients[[2]]
        }
      })
    MCAv[m, , ] <- apply(Y, c(2, 3),
      function(time_series) {
        if (any(is.na(time_series))) {
          NA
        } else {
          lm(time_series ~ ECv_std[m, ])$coefficients[[2]]
        }
      })
  }

  invisible(list(SC = SV2[1:nm], SCFs = SCF, RUVs = RUV, ECs_U = ECu_std, 
                 MCAs_U = MCAu, ECs_V = ECv_std, MCAs_V = MCAv)) 
  #               CORRs_U = CORRU, CORRs_V = CORRV))
}
