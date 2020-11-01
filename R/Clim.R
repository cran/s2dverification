#'Computes Bias Corrected Climatologies
#'
#'This function computes only per-pair climatologies from the experimental 
#'and observational matrices output from \code{Load()}.
#'To compute plain climatologies from only experimental or observational 
#'data from \code{Load()}, the following code can be used:\cr 
#'\code{clim <- array(apply(obs_data, c(1, 4, 5, 6), mean),}\cr
#'\code{              dim = dim(obs_datta)[-c(2, 3)])}\cr
#'The function \code{Clim()} computes per-pair climatologies using one of the 
#'following methods:
#'\enumerate{
#'  \item{per-pair method (Garcia-Serrano and Doblas-Reyes, CD, 2012)}
#'  \item{Kharin method (Karin et al, GRL, 2012)}
#'  \item{Fuckar method (Fuckar et al, GRL, 2014)}
#'}
#'\code{Clim()} computes climatologies using the startdates covered by the 
#'whole experiments/observational data sets. The startdates not available for 
#'all the data (model and obs) are excluded when computing the climatologies.
#'
#'@param var_exp Model data: c(nmod/nexp, nmemb/nparam, nsdates, nltime) up to 
#'  c(nmod/nexp, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon).
#'@param var_obs Observational data: c(nobs, nmemb, nsdates, nltime) up to 
#'  c(nobs, nmemb, nsdates, nltime, nlevel, nlat, nlon).
#'@param memb TRUE/FALSE (1 climatology for each member). Default = TRUE.
#'@param kharin TRUE/FALSE (if Kharin method is applied or not). 
#'  Default = FALSE.
#'@param NDV TRUE/FALSE (if Fuckar method is applied or not). Default = FALSE.
#'
#'@return
#'\item{clim_exp}{Array with same dimensions as var_exp except the third 
#'  (starting dates) and, depending on the parameters, the second (members), 
#'  which disappear.}
#'\item{clim_obs}{Array with same dimensions as var_obs except the third 
#'  (starting dates) and, depending on the parameters, the second (members), 
#'  which disappear.}
#'@keywords datagen
#'@author History:\cr
#'  0.9 - 2011-03 (V. Guemas, \email{virginie.guemas@@ic3.cat}) - Original code\cr
#'  1.0 - 2013-09 (N. Manubens, \email{nicolau.manubens@@ic3.cat}) - Formatting to R CRAN
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
#'clim <- Clim(sampleData$mod, sampleData$obs)
#'  \donttest{
#'PlotClim(clim$clim_exp, clim$clim_obs, 
#'         toptitle = paste('sea surface temperature climatologies'), 
#'         ytitle = 'K', monini = 11, listexp = c('CMIP5 IC3'), 
#'         listobs = c('ERSST'), biglab = FALSE, fileout = 'tos_clim.eps')
#'  }
#'@export
Clim <- function(var_exp, var_obs, memb = TRUE, kharin = FALSE, NDV = FALSE) {
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
  #  Find common points to compute climatologies 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #

  na_array <- array(0, dim = dim(var_exp)[3:7])

  for (i_dat in 1:dimexp[1]) {
    for (i_memb in 1:dimexp[2]) {
      na_array <- na_array + as.numeric(is.na(var_exp[i_dat, i_memb, , , , , ]))
    }
  }
  for (i_dat in 1:dimobs[1]) {
    for (i_memb in 1:dimobs[2]) {
      na_array <- na_array + as.numeric(is.na(var_obs[i_dat, i_memb, , , , , ]))
    }
  }
  na_array <- as.logical(na_array)  # TRUE is NA, FALSE is not

  for (i_dat in 1:dimexp[1]) {
    for (i_memb in 1:dimexp[2]) {      
      asd <- var_exp[i_dat, i_memb, , , , , ]
      asd[which(na_array)] <- NA
      var_exp[i_dat, i_memb, , , , , ] <- asd
    }
  }
  for (i_dat in 1:dimobs[1]) {
    for (i_memb in 1:dimobs[2]) {
      asd <- var_obs[i_dat, i_memb, , , , , ]
      asd[which(na_array)] <- NA
      var_obs[i_dat, i_memb, , , , , ] <- asd
    }
  }

#  nan <- MeanListDim(var_exp, dims = c(1:2, 5:7), narm = FALSE) + 
#         MeanListDim(var_obs, dims = c(1:2, 5:7), narm = FALSE)
#
#  for (jdate in 1:dimexp[3]) {
#    for (jtime in 1:dimexp[4]) {
#      if (is.na(nan[jdate, jtime])) {
#        var_exp[, , jdate, jtime, , , ] <- NA
#        var_obs[, , jdate, jtime, , , ] <- NA
#      }
#    }
#  }



  #
  #  Compute climatologies 
  # ~~~~~~~~~~~~~~~~~~~~~~~
  #
  out_clim_obs <- Mean1Dim(var_obs, posdim = 3, narm = TRUE)
  dim_clim_obs <- dimobs[-3]

  if (kharin == TRUE) {
    tmp_obs <- Trend(var_obs, posTR = 3)
    tmp_exp <- Trend(var_exp, posTR = 3)
    trend_obs <- array(dim = dim(var_exp))
    trend_exp <- array(dim = dim(var_exp))
    for (jdate in 1:dimexp[3]) {
      trend_exp[, , jdate, , , , ] <- tmp_exp$trend[, , 4, , , , ] +
                                       jdate * tmp_exp$trend[, , 2, , , , ]
      tmp_obs2 <- MeanListDim(tmp_obs$trend,c(2, 1))
      trend_obs[, , jdate, , , , ] <- InsertDim(InsertDim(tmp_obs2[4, , , , ] + 
                                      jdate * tmp_obs2[2, , , , ], 1, dimexp[1]),
                                      2, dimexp[2])
    }
    out_clim_exp <- trend_exp - trend_obs + InsertDim(InsertDim(InsertDim(
                    MeanListDim(out_clim_obs, c(2, 1)), 1, dimexp[1]), 2, 
                    dimexp[2]), 3, dimexp[3])
    dim_clim_exp <- dimexp
  } else if (NDV == TRUE) {
    iniobs <- InsertDim(SelIndices(var_obs, 4, c(1, 1)), 4, dimobs[4])
    iniexp <- InsertDim(SelIndices(var_exp, 4, c(1, 1)), 4, dimexp[4])
    tmp_obs <- Regression(var_obs, iniobs, posREG = 3)
    tmp_exp <- Regression(var_exp, iniexp, posREG = 3)
    reg_obs <- array(dim = dim(var_exp))
    reg_exp <- array(dim = dim(var_exp))
    for (jdate in 1:dimexp[3]) {
      reg_exp[, , jdate, , , , ] <- tmp_exp$regression[, , 4, , , , ] +
                                     iniexp[, , jdate, , , , ] * 
                                    tmp_exp$regression[, , 2, , , , ]
      tmp_obs2 <- MeanListDim(tmp_obs$regression,c(2, 1))
      reg_obs[, , jdate, , , , ] <- InsertDim(InsertDim(tmp_obs2[4, , , , ] +
                                     MeanListDim(iniobs,c(2, 1))[jdate, , , , ] *
                                     tmp_obs2[2, , , , ], 1, dimexp[1]), 
                                    2, dimexp[2])
    }
    out_clim_exp <- reg_exp - reg_obs + InsertDim(InsertDim(InsertDim(
                    MeanListDim(out_clim_obs, c(2, 1)), 1, dimexp[1]), 2, 
                    dimexp[2]), 3, dimexp[3])
    dim_clim_exp <- dimexp
  } else {
    out_clim_exp <- Mean1Dim(var_exp, posdim = 3, narm = TRUE)
    dim_clim_exp <- dimexp[-3]
  }

  if (memb != TRUE) {
    out_clim_obs <- Mean1Dim(out_clim_obs, posdim = 2, narm = TRUE) 
    out_clim_exp <- Mean1Dim(out_clim_exp, posdim = 2, narm = TRUE)
    dim_clim_exp <- dim_clim_exp[-2]
    dim_clim_obs <- dim_clim_obs[-2]
  }

  #
  #  Reduce the number of dimensions to the original one 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 

  clim_exp <- array(dim = dim_clim_exp)
  clim_exp[] <- out_clim_exp
  clim_obs <- array(dim = dim_clim_obs)
  clim_obs[] <- out_clim_obs

  #
  #  Outputs
  # ~~~~~~~~~
  #
  invisible(list(clim_exp = clim_exp, clim_obs = clim_obs))
}
