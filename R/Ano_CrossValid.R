#'Computes Anomalies In Cross-Validation Mode
#'
#'Computes the anomalies from the arrays of the experimental and observational 
#'data output from \code{load()} by subtracting the climatologies computed 
#'with a cross-validation technique and a per-pair method.
#'
#'@param var_exp Model data:\cr
#'c(nmod/nexp, nmemb/nparam, nsdates, nltime) up to \cr
#'  c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon).
#'@param var_obs Observational data: \cr
#'c(nobs, nmemb, nsdates, nltime) up to \cr
#'  c(nobs, nmemb, nsdates, nltime, nlevel, nlat, nlon).
#'@param memb TRUE/FALSE (1 climatology for each member/1 climatology 
#'  averaging all the members). Default = TRUE.
#'
#'@return
#'  \item{$ano_exp}{Matrix with same dimensions as var_exp}
#'  \item{$ano_obs}{Matrix with same dimensions as var_obs}
#'
#'@keywords datagen
#'@author History:\cr
#'  0.1 - 2011-12 (V. Guemas) - Original code\cr
#'  1.0 - 2013-09 (N. Manubens) - Formatting to CRAN
#'
#'@examples 
#'# Load sample data as in Load() example:
#'example(Load)
#'anomalies <- Ano_CrossValid(sampleData$mod, sampleData$obs)
#'  \donttest{
#'PlotAno(anomalies$ano_exp, anomalies$ano_obs, startDates, 
#'        toptitle = paste('anomalies'), ytitle = c('K', 'K', 'K'), 
#'        legends = 'ERSST', biglab = FALSE, fileout = 'tos_ano_crossvalid.eps')
#'  }
#'@export
Ano_CrossValid <- function(var_exp, var_obs, memb = TRUE) {
  #
  #  Enlarge the number of dimensions of var_exp and var_obs to 7 if necessary
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimexp <- dim(var_exp)
  dimobs <- dim(var_obs)
  if (length(dimexp) < 4 | length(dimobs) < 4) {
    stop("At least 4 dim needed : c(nexp/nobs, nmemb, nsdates, nltime)")
  }
  for (jn in 3:max(length(dimexp), length(dimobs))) {
    if (dimexp[jn] != dimobs[jn]) { 
      stop("Wrong input dimensions")
    }
  }
  var_exp <- Enlarge(var_exp, 7)
  var_obs <- Enlarge(var_obs, 7)
  #
  #  Computation of anomalies in cross-validation mode
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  enl_ano_exp <- array(dim = dim(var_exp))
  enl_ano_obs <- array(dim = dim(var_obs))
  for (jsdat in 1:dimexp[3]) {
    tmp1 <- array(var_exp[, , jsdat, , , , ], dim = dim(var_exp)[-3])
    tmp2 <- array(var_obs[, , jsdat, , , , ], dim = dim(var_obs)[-3])
    tmp3 <- array(var_exp[, , -jsdat, , , , ], dim = c(dimexp[1:2], 
                                                       dimexp[3] - 1, 
                                                       dim(var_exp)[4:7]))
    tmp4 <- array(var_obs[, , -jsdat, , , , ], dim = c(dimobs[1:2], 
                                                       dimobs[3] - 1, 
                                                       dim(var_obs)[4:7]))
    tmp <- Clim(tmp3, tmp4, memb)
    if (memb) { 
      clim_exp <- tmp$clim_exp
      clim_obs <- tmp$clim_obs
    } else {
      clim_exp <- InsertDim(tmp$clim_exp, 2, dimexp[2]) 
      clim_obs <- InsertDim(tmp$clim_obs, 2, dimobs[2]) 
    }
    enl_ano_exp[, , jsdat, , , , ] <- tmp1 - clim_exp
    enl_ano_obs[, , jsdat, , , , ] <- tmp2 - clim_obs  
  } 
  #
  #  Reduce the number of dimensions to the original one 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  ano_exp <- array(dim = dimexp)
  ano_obs <- array(dim = dimobs)
  ano_exp[] <- enl_ano_exp
  ano_obs[] <- enl_ano_obs
  #
  #  Outputs
  # ~~~~~~~~~
  #
  invisible(list(ano_exp = ano_exp, ano_obs = ano_obs))
}
