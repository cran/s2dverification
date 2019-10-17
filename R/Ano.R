#'Computes Forecast or Observed Anomalies
#'
#'This function computes anomalies from any experimental or observational 
#'matrix output from \code{Load()} and their climatologies output from 
#'\code{Clim()}.
#'
#'@param var Model or observational data:\cr 
#'  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime) up to \cr
#'  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)\cr
#'@param clim Climatologies from clim: c(nmod/nexp/nobs, nltime) \cr
#'  up to c(nmod/nexp/nobs, nltime, nlevel, nlat, nlon) or \cr
#'  c(nmod/nexp/nobs, nmemb/nparam, nltime) up to \cr
#'  c(nmod/nexp/nobs, nmemb/nparam, nltime, nlevel, nlat, nlon) or \cr
#'  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime) up to \cr
#'  c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon) \cr
#'  depending on the options provided to \code{Clim()}.
#'
#'@return Array with same dimensions as 'var'.
#'@keywords datagen
#'@author History:\cr
#'  0.1 - 2012-03 (V. Guemas, \email{virginie.guemas@@ic3.cat}) - Original code\cr
#'  1.0 - 2013-09 (N. Manubens, \email{nicolau.manubens@@ic3.cat}) - Formatting to R CRAN
#'
#'@examples 
#'# Load sample data as in Load() example:
#'example(Load)
#'clim <- Clim(sampleData$mod, sampleData$obs)
#'ano_exp <- Ano(sampleData$mod, clim$clim_exp)
#'ano_obs <- Ano(sampleData$obs, clim$clim_obs)
#'runmean_nb_months <- 12
#'dim_to_smooth <- 4  # Smooth along lead-times
#'smooth_ano_exp <- Smoothing(ano_exp, runmean_nb_months, dim_to_smooth)
#'smooth_ano_obs <- Smoothing(ano_obs, runmean_nb_months, dim_to_smooth)
#'  \donttest{
#'PlotAno(smooth_ano_exp, smooth_ano_obs, startDates, 
#'        toptitle = paste('smoothed anomalies'), ytitle = c('K', 'K', 'K'), 
#'        legends = 'ERSST', biglab = FALSE, fileout = 'tos_ano.eps')
#'  }
#'@export
Ano <- function(var, clim) {
  #
  #  Duplicate clim dimensions to heve same dimensions as var
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimvar <- dim(var)

  if (length(dim(clim)) <= 2) {
    clim <- InsertDim(clim, 2, dimvar[2])
  }

  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[2] != dimvar[2])) {
    clim <- InsertDim(clim, 2, dimvar[2])
  }
  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[2] == dimvar[2])) {
    if (is.null(names(dim(clim))) | is.null(names(dimvar))) {
      stop('Provide dimension names on parameter \'var\' and \'clim\' to avoid ambiguity.')
    } else {
      if (names(dim(clim)[2]) != names(dimvar[2])) {
        clim <- InsertDim(clim, 2, dimvar[2])
      }
    }
  }

  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[3] != dimvar[3])) {
    clim <- InsertDim(clim, 3, dimvar[3])
  }
  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[3] == dimvar[3])) {
    if (is.null(names(dim(clim))) | is.null(names(dimvar))) {
      stop('Provide dimension names on parameter \'var\' and \'clim\' to avoid ambiguity.')
    } else {
      if (names(dim(clim)[3]) != names(dimvar[3])) {
        clim <- InsertDim(clim, 3, dimvar[3])
      }
    }
  }

  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[4] != dimvar[4])) {
    clim <- InsertDim(clim, 4, dimvar[4])
  }
  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[4] == dimvar[4])) {
    if (is.null(names(dim(clim))) | is.null(names(dimvar))) {
      stop('Provide dimension names on parameter \'var\' and \'clim\' to avoid ambiguity.')
    } else {
      if (names(dim(clim)[4]) != names(dimvar[4])) {
        clim <- InsertDim(clim, 4, dimvar[4])
      }
    }
  }

  #
  #  Raw anomalies
  # ~~~~~~~~~~~~~~~
  #
  ano <- var - clim

  #
  #  Outputs
  # ~~~~~~~~~
  #
  ano
}
