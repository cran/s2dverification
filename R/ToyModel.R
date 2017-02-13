ToyModel <- function(alpha = 0.1, beta = 0.4, gamma = 1, sig = 1, 
  trend = 0, nstartd = 30, nleadt = 4, nmemb = 10, obsini = NULL, 
  fxerr = NULL) {
  # Check alpha, beta, gamma, sig, trend, nstartd, nleadt, nmemb
  if (any(!is.numeric(c(alpha, beta, gamma, sig, trend, nstartd, 
                      nleadt, nmemb)))) {
    stop(paste("Parameters 'alpha', 'beta', 'gamma', 'sig', 'trend', 'nstartd',",
               "nleadt and nmemb must be numeric."))
  }
  nstartd <- round(nstartd)
  nleadt <- round(nleadt)
  nmemb <- round(nmemb)

  # Check obsini
  if (!is.null(obsini)) {
    if (!is.numeric(obsini) || !is.array(obsini)) {
      stop("Parameter 'obsini' must be a numeric array.")
    }
    if (length(dim(obsini)) != 4) {
      stop("Parameter 'obsini' must be an array with dimensions c(1, 1, nleadt, nstartd).")
    }
    if (dim(obsini)[3] != nstartd || dim(obsini)[4] != nleadt) {
      stop(paste0("The dimensions of parameter 'obsini' and the parameters 'nleadt' and 'nstartd' must match:\n  dim(obsini) = c(", 
           dim(obsini)[3], ", ", dim(obsini)[4], ")\n  nstartd = ", nstartd, "  nleadt = ", nleadt))
    }
  }

  # Check fxerr
  if (!is.null(fxerr)) {
    if (!is.numeric(fxerr)) {
      stop("Parameter 'fxerr' must be numeric.")
    }
  }

  time <- seq(1, nstartd)  # time vector, generated from 1 -> nstard
  
  if (!(sig^2 - alpha^2 - beta^2 > 0)) {
    stop("Model variability not constrained: respect condition \"sig^2-alpha^2-beta^2 > 0\")")
  }
  
  if (nstartd < 0) {
    stop("Number of start dates must be positive")
  }
  
  if (nleadt < 0) {
    stop("Number of lead-times must be positive")
  }
  
  if (nmemb < 0) {
    stop("Number of members must be positive")
  }
  
  
  if (!is.null(obsini)) {
    # Observations provided by input compute forecast parameters
    # from observations
    obs_ano <- obsini
  } else {
    # Create observations for selected period
    obs_ano <- array(rnorm(nleadt * nstartd, mean = 0, sd = sig), 
      dim = c(1, 1, nstartd, nleadt))
    obs_trend <- array(t(time) * rep(trend, times = nstartd), 
      , dim = c(1, 1, nstartd, nleadt))
    obs <- obs_ano + obs_trend
    
    trend <- rep(c(trend), times = nleadt)  # same trend in all lead-times
    sig <- rep(c(sig), times = nleadt)  # same variability in all lead-times
  }
  
  forecast <- array(dim = c(length(gamma), nmemb, nstartd, 
    nleadt))
  
  # Allocate observations and forecast according to
  # s2dverification standards
  for (j in 1:nstartd) {
    for (f in 1:nleadt) {
      for (g in 1:length(gamma)) {
        # Generate forecasts with different trends
        auto_term <- alpha * obs_ano[1, 1, j, f]  # Predictability term
        if (is.numeric(fxerr)) {
          conf_term <- fxerr  # Forecast error term, fixed by input
        } else {
          conf_term <- rnorm(1, mean = 0, sd = beta)  # Forecast error term, random
        }
        trend_term <- gamma[g] * trend[f] * j  # Trend term
        var_corr <- rnorm(nmemb, mean = 0, sd = sqrt(sig[f] - 
          alpha^2 - beta^2))
        forecast[g, , j, f] <- matrix(auto_term, c(nmemb, 
          1)) + matrix(conf_term, c(nmemb, 1)) + matrix(trend_term, 
          c(nmemb, 1)) + var_corr
      }
    }
  }
  
  list(mod = forecast, obs = obs_ano)
}
