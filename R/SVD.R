#'Single Value Decomposition (Maximum Covariance Analysis)
#'
#'Computes a Maximum Covariance Analysis (MCA) between vary and varx, both
#'of dimensions c(n. of time steps, n. of latitudes, n. of longitudes), each
#'over a region of interest, e.g.: prlr over Europe and tos over North Atlantic.
#'The input fields are latitude-weighted by default (can be adjustable via 
#'\code{weight}).\cr
#'Returns a vector of squared covariance fraction (SCFs) explained by 
#'each pair of covariability modes, a vector of correlation coefficient 
#'(RUVs) between expansion coefficients (ECs) that measures their linear 
#'relationship, and a set of regression (MCAs) associated with the 
#'covariability modes (ECs). Note that MCAs are 'homogeneous' patterns obtained 
#'as regression/correlation between each field (predictor, predictand) 
#'and its expansion coefficient.\cr
#'The MCA is computed by default with the covariance matrix. It can be computed
#'with the correlation matrix by setting \code{corr = TRUE}.
#'
#'@param vary Array containing the anomalies field for the predictor. The 
#'  expected dimensions are c(n. of time steps, n. of latitudes, n. of 
#'  longitudes).
#'@param varx Array containing the anomalies field for the predictand. The 
#'  expected dimensions are c(n. of time steps, n. of latitudes, n. of 
#'  longitudes).
#'@param laty Vector of latitudes of the array \code{vary}. Only required if 
#'  \code{weight = TRUE}.
#'@param latx Vector of latitudes of the array \code{varx}. Only required if 
#'  \code{weight = TRUE}.
#'@param nmodes Number of ECs/MCAs/modes retained and provided in the outputs.
#'@param corr Whether to compute the MCA over a covariance matrix (FALSE) or 
#'  a correlation matrix (TRUE).
#'@param weight Whether to apply latitude weights on the input fields or not. 
#'  TRUE by default.
#'
#'@return 
#'  \item{$SC}{
#'Vector of squared covariance (n. of modes).
#'  }
#'  \item{$SCFs}{
#'    Vector of squared covariance fractions (n. of modes).
#'  }
#'  \item{$RUVs}{
#'    Vector of correlations between expansion coefficients (n. of modes).
#'  }
#'  \item{$ECs_U}{
#'    Array of expansion coefficients of predictor field (n. of time steps, 
#'    n. of modes).
#'  }
#'  \item{$MCAs_U}{
#'    Array of covariability patterns of predictor field (c(dim), n. of modes).
#'  }
#'  \item{$ECs_V}{
#'    Array of expansion coefficients of predictand field (n. of time steps, 
#'    n. of modes).
#'  }
#'  \item{$MCAs_V}{
#'    Array of covariability patterns of predictand field (c(dim), n. of modes).
#'  }
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2010-09  (J.-G. Serrano)  -  Original code\cr
#'1.0  -  2016-04  (N. Manubens)  -  Formatting to R CRAN
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
#'# This example computes the ECs and MCAs along forecast horizons and plots the 
#'# one that explains the greatest amount of variability. The example data is 
#'# very low resolution so it does not make a lot of sense.
#'ano <- Ano_CrossValid(sampleData$mod, sampleData$obs)
#'mca <- SVD(Mean1Dim(ano$ano_exp, 2)[1, , 1, , ], 
#'           Mean1Dim(ano$ano_obs, 2)[1, , 1, , ], 
#'           sampleData$lat, sampleData$lat)
#'PlotEquiMap(mca$MCAs_U[1, , ], sampleData$lon, sampleData$lat)
#'plot(mca$ECs_U[1, ])
#'PlotEquiMap(mca$MCAs_V[1, , ], sampleData$lon, sampleData$lat)
#'plot(mca$ECs_V[1, ])
#'
#'@importFrom stats sd cor lm
#'@export
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
