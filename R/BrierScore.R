#'Compute Brier Score And Its Decomposition And Brier Skill Score
#'
#'Computes the Brier score (BS) and the components of its standard 
#'decomposition as well with the two within-bin components described in 
#'Stephenson et al., (2008). It also returns the bias-corrected decomposition 
#'of the BS (Ferro and Fricker, 2012). BSS having the climatology as the 
#'reference forecast. \cr\cr
#'.BrierScore provides the same functionality, but taking a matrix of ensemble 
#'members (exp) as input.
#'
#'@param obs Vector of binary observations (1 or 0).
#'@param pred Vector of probablistic predictions with values in the range [0,1].
#'@param thresholds Values used to bin the forecasts. By default the bins are 
#'  {[0,0.1), [0.1, 0.2), ... [0.9, 1]}.
#'@param exp Matrix of predictions with values in the range [0,1] for the 
#'  .BrierScore function
#'
#'@return Both BrierScore and .Brier score provide the same outputs:
#'\itemize{ 
#'  \item{$rel}{standard reliability}
#'  \item{$res}{standard resolution}
#'  \item{$unc}{standard uncertainty}  
#'  \item{$bs}{Brier score}
#'  \item{$bs_check_res}{rel-res+unc}
#'  \item{$bss_res}{res-rel/unc}
#'  \item{$gres}{generalized resolution}
#'  \item{$bs_check_gres}{rel-gres+unc}
#'  \item{$bss_gres}{gres-rel/unc}
#'  \item{$rel_bias_corrected}{bias-corrected rel}
#'  \item{$gres_bias_corrected}{bias-corrected gres}
#'  \item{$unc_bias_corrected}{bias-corrected unc}
#'  \item{$bss_bias_corrected}{gres_bias_corrected-rel_bias_corrected/unc_bias_corrected}
#'  \item{$nk}{number of forecast in each bin}
#'  \item{$fkbar}{average probability of each bin}
#'  \item{$okbar}{relative frequency that the observed event occurred}
#'  \item{$bins}{bins used}
#'  \item{$pred}{values with which the forecasts are verified}
#'  \item{$obs}{probability forecasts of the event}
#'}
#'
#'@references
#'Wilks (2006) Statistical Methods in the Atmospheric Sciences.\cr
#'Stephenson et al. (2008). Two extra components in the Brier score decomposition. 
#'  Weather and Forecasting, 23: 752-757.\cr
#'Ferro and Fricker (2012). A bias-corrected decomposition of the BS. 
#'  Quarterly Journal of the Royal Meteorological Society, DOI: 10.1002/qj.1924.
#'
#'@keywords datagen
#'@author History:\cr
#'  0.1 - 2012-04 (L. Rodrigues, \email{lrodrigues@@ic3.cat}) - Original code\cr
#'  0.2 - 2017-02 (A. Hunter, \email{alasdair.hunter@@bsc.es}) - Adapted to veriApply()
#'
#'@examples
#'# Minimalist examples with BrierScore
#'a <- runif(10)
#'b <- round(a)
#'x <- BrierScore(b, a)
#'x$bs - x$bs_check_res
#'x$bs - x$bs_check_gres
#'x$rel_bias_corrected - x$gres_bias_corrected + x$unc_bias_corrected
#'  \dontrun{
#'a <- runif(10)
#'b <- cbind(round(a),round(a)) # matrix containing 2 identical ensemble members...
#'x2 <- BrierScore(a, b)
#'  }
#'
#'# Example of BrierScore using UltimateBrier
#'# See ?UltimateBrier for more information
#'example(Load)
#'clim <- Clim(sampleData$mod, sampleData$obs)
#'ano_exp <- Ano(sampleData$mod, clim$clim_exp)
#'ano_obs <- Ano(sampleData$obs, clim$clim_obs)
#'bs <- UltimateBrier(ano_exp, ano_obs, thr = c(1/3, 2/3))
#'
#'  \dontrun{
#'# Example of .BrierScore with veriApply
#'require(easyVerification)
#'BrierScore2 <- s2dverification:::.BrierScore
#'bins_ano_exp <- ProbBins(ano_exp, thr = c(1/3, 2/3), posdates = 3, posdim = 2)
#'bins_ano_obs <- ProbBins(ano_obs, thr = c(1/3, 2/3), posdates = 3, posdim = 2)
#'bs2 <- veriApply("BrierScore2", bins_ano_exp, Mean1Dim(bins_ano_ob,s 3), 
#'                 tdim = 2, ensdim = 3)
#'  }
#'@rdname BrierScore
#'@export
BrierScore <- function(obs, pred, thresholds = seq(0, 1, 0.1)) {
  if (max(pred) > 1 | min(pred) < 0) {
    stop("Predictions outside [0,1] range. Are you certain this is a probability forecast? \n")
  } else if (max(obs) != 1 & min(obs) != 0) {
    .message("Binary events must be either 0 or 1. Are you certain this is a binary event? ")
  } else {
    nbins <- length(thresholds) - 1  # Number of bins
    n <- length(pred)
    bins <- as.list(paste("bin", 1:nbins,sep = ""))
    for (i in 1:nbins) {
      if (i == nbins) {
        bins[[i]] <- list(which(pred >= thresholds[i] & pred <= thresholds[i + 1]))
      } else {
        bins[[i]] <- list(which(pred >= thresholds[i] & pred < thresholds[i + 1]))
      }
    }
    
    fkbar <- okbar <- nk <- array(0, dim = nbins)
    for (i in 1:nbins) {
      nk[i] <- length(bins[[i]][[1]])
      fkbar[i] <- sum(pred[bins[[i]][[1]]]) / nk[i]
      okbar[i] <- sum(obs[bins[[i]][[1]]]) / nk[i]
    }
    
    obar <- sum(obs) / length(obs)
    relsum <- ressum <- term1 <- term2 <- 0
    for (i in 1:nbins) {
      if (nk[i] > 0) {
        relsum <- relsum + nk[i] * (fkbar[i] - okbar[i])^2
        ressum <- ressum + nk[i] * (okbar[i] - obar)^2
        for (j in 1:nk[i]) {
          term1 <- term1 + (pred[bins[[i]][[1]][j]] - fkbar[i])^2
          term2 <- term2 + (pred[bins[[i]][[1]][j]] - fkbar[i]) * (obs[bins[[i]][[1]][j]] - okbar[i])
        }
      }
    }
    rel <- relsum / n
    res <- ressum / n
    unc <- obar * (1 - obar)
    bs <- sum((pred - obs)^2) / n
    bs_check_res <- rel - res + unc
    bss_res <- (res - rel) / unc
    gres <- res - term1 * (1 / n) + term2 * (2 / n)   # Generalized resolution
    bs_check_gres <- rel - gres + unc                 # BS using GRES
    bss_gres <- (gres - rel) / unc                    # BSS using GRES
    
    #
    # Estimating the bias-corrected components of the BS 
    #
    term3 <- array(0, nbins)
    for (i in 1:nbins) {
      term3[i] <- (nk[i] / (nk[i] - 1)) * okbar[i] * (1 - okbar[i])
    }
    term_a <- sum(term3,  na.rm = T) / n
    term_b <- (obar * (1 - obar)) / (n - 1)
    rel_bias_corrected <- rel - term_a
    gres_bias_corrected <- gres - term_a + term_b
    if (rel_bias_corrected < 0 || gres_bias_corrected < 0) {
      rel_bias_corrected2 <- max(rel_bias_corrected, rel_bias_corrected - gres_bias_corrected, 0)
      gres_bias_corrected2 <- max(gres_bias_corrected, gres_bias_corrected - rel_bias_corrected, 0)
      rel_bias_corrected <- rel_bias_corrected2
      gres_bias_corrected <- gres_bias_corrected2
    }
    unc_bias_corrected <- unc + term_b
    bss_bias_corrected <- (gres_bias_corrected - rel_bias_corrected) / unc_bias_corrected
    
    #if (round(bs, 8) == round(bs_check_gres, 8) & round(bs_check_gres, 8) == round((rel_bias_corrected - gres_bias_corrected + unc_bias_corrected), 8)) {
    #  cat("No error found  \ n")
    #  cat("BS = REL - GRES + UNC = REL_lessbias - GRES_lessbias + UNC_lessbias  \ n")
    #}
    
    invisible(list(rel = rel, res = res, unc = unc, bs = bs, bs_check_res = bs_check_res, bss_res = bss_res, gres = gres, bs_check_gres = bs_check_gres, bss_gres = bss_gres, rel_bias_corrected = rel_bias_corrected, gres_bias_corrected = gres_bias_corrected, unc_bias_corrected = unc_bias_corrected, bss_bias_corrected = bss_bias_corrected, nk = nk, fkbar = fkbar, okbar = okbar, bins = bins, pred = pred, obs = obs))
  }
}

#'@rdname BrierScore
#'@export
.BrierScore <- function(exp, obs, thresholds = seq(0, 1, 0.1)) {
  if (max(exp) > 1 || min(exp) < 0) {
    stop("Parameter 'exp' contains predictions outside [0,1] range. Are you certain this is a probability forecast?")
  } else if (max(obs) != 1 && min(obs) != 0) {
    .message("Binary events in 'obs' must be either 0 or 1. Are you certain this is a binary event?")
  } else {
    nbins <- length(thresholds) - 1  # Number of bins
    n <- dim(exp)[1]                 # Number of observations
    ens_mean <- rowMeans(exp, na.rm = TRUE)
    n.ens <- seq(1, dim(exp)[2], 1)  # Number of ensemble members
    bins <- as.list(paste("bin", 1:nbins, sep = ""))
    for (i in 1:nbins) {
      if (i == nbins) {
        bins[[i]] <- list(which(ens_mean >= thresholds[i] & ens_mean <= thresholds[i + 1])) 
      } else {
        bins[[i]] <- list(which(ens_mean >= thresholds[i] & ens_mean < thresholds[i + 1])) 
      }
    }

    fkbar <- okbar <- nk <- array(0, dim = nbins)
    for (i in 1:nbins) {
      nk[i] <- length(bins[[i]][[1]])
      fkbar[i] <- sum(ens_mean[bins[[i]][[1]]]) / nk[i]
      okbar[i] <- sum(obs[bins[[i]][[1]]]) / nk[i]
    }

    fkbar[fkbar == Inf] <- 0 
    okbar[is.nan(okbar)] <- 0
    obar <- sum(obs) / length(obs)
    relsum <- ressum <- relsum1 <- ressum1 <- term1 <- term1a <- term2 <- term2a <- 0

    for (i in 1:nbins) {
      if (nk[i] > 0) {
        relsum <- relsum + nk[i] * (fkbar[i] - okbar[i]) ^ 2
        ressum <- ressum + nk[i] * (okbar[i] - obar) ^ 2
        for (j in 1:nk[i]) {
          term1 <- term1 + (ens_mean[bins[[i]][[1]][j]] - fkbar[i]) ^ 2
          term2 <- term2 + (ens_mean[bins[[i]][[1]][j]] - fkbar[i]) * (obs[bins[[i]][[1]][j]] - okbar[i])
        }
      }
    }
  }
  rel <- relsum / n
  res <- ressum / n
  unc <- obar * (1 - obar)
  #bs <- apply(ens, MARGIN = 2, FUN = function(x) sum((x - obs)^2) / n)
  bs <- sum((rowMeans(exp, na.rm = T) - obs) ^ 2) / n
  bs_check_res <- rel - res + unc
  bss_res <- (res - rel) / unc
  gres <- res - term1 * (1 / n) + term2 * (2 / n)   # Generalized resolution
  bs_check_gres <- rel - gres + unc                 # BS using GRES
  bss_gres <- (gres - rel) / unc                    # BSS using GRES
  
  #
  # Estimating the bias-corrected components of the BS 
  #
  term3 <- array(0, nbins)
  for (i in 1:nbins) {
    term3[i] <- (nk[i] / (nk[i] - 1)) * okbar[i] * (1 - okbar[i])
  }
  term_a <- sum(term3,  na.rm = T) / n
  term_b <- (obar * (1 - obar)) / (n - 1)
  rel_bias_corrected <- rel - term_a
  gres_bias_corrected <- gres - term_a + term_b
  if (rel_bias_corrected < 0 || gres_bias_corrected < 0) {
    rel_bias_corrected2 <- max(rel_bias_corrected, rel_bias_corrected - gres_bias_corrected, 0)
    gres_bias_corrected2 <- max(gres_bias_corrected, gres_bias_corrected - rel_bias_corrected, 0)
    rel_bias_corrected <- rel_bias_corrected2
    gres_bias_corrected <- gres_bias_corrected2
  }
  unc_bias_corrected <- unc + term_b
  bss_bias_corrected <- (gres_bias_corrected - rel_bias_corrected) / unc_bias_corrected
    
  #if (round(bs, 8) == round(bs_check_gres, 8) & round(bs_check_gres, 8) == round((rel_bias_corrected - gres_bias_corrected + unc_bias_corrected), 8)) {
  #  cat("No error found  \ n")
  #  cat("BS = REL - GRES + UNC = REL_lessbias - GRES_lessbias + UNC_lessbias  \ n")
  #}
    
  invisible(list(rel = rel, res = res, unc = unc, bs = bs, 
                 bs_check_res = bs_check_res, bss_res = bss_res, gres = gres, 
                 bs_check_gres = bs_check_gres, bss_gres = bss_gres, 
                 rel_bias_corrected = rel_bias_corrected, 
                 gres_bias_corrected = gres_bias_corrected, 
                 unc_bias_corrected = unc_bias_corrected, 
                 bss_bias_corrected = bss_bias_corrected))
}
