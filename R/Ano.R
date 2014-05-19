Ano <- function(var, clim) {
  # This function computes raw anomalies from experimental or observational 
  # matrix output from load() and their climatologies output from clim()
  #
  # Args:
  #   var: Model or observational data:
  #        c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime) up to
  #        c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime, nlevel, nlat, nlon)
  #   clim: Climatologies from clim:
  #         c(nmod/nexp/nobs, nmemb/nparam, nltime) up to
  #         c(nmod/nexp/nobs, nmemb/nparam, nltime, nlevel, nlat, nlon)
  #         or
  #         c(nmod/nexp/nobs, nltime) up to
  #         c(nmod/nexp/nobs, nltime, nlevel, nlat, nlon)
  #
  # Returns:
  #   Matrix with same dimensions as var
  #
  # History:
  #   1.0  #  2011-03  (V. Guemas, vguemas@ic3.cat)  #  Original code
  
  #
  #  Duplicate clim dimensions to heve same dimensions as var
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimvar <- dim(var)
  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[2] != dimvar[2])) {
    clim <- InsertDim(clim, 2, dimvar[2])
  }
  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[3] != dimvar[3])) {
    clim <- InsertDim(clim, 3, dimvar[3])
  }
  if ((length(dimvar) > length(dim(clim))) & (dim(clim)[4] != dimvar[4])) {
    clim <- InsertDim(clim, 4, dimvar[4])
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
