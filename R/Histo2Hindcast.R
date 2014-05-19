Histo2Hindcast <- function(varin, sdatesin, sdatesout, nleadtimesout) {
  # This function reorganizes a long run (historical typically) with only one 
  # start date into chunks corresponding to a set of start dates. The expected 
  # input structure is the one output from load with 4 to 7 dimensions.
  #
  # Args:
  #   varin: Input model or observational data:
  #            c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime) up to
  #            c(nmod/nexp/nobs, nmemb/nparam, nsdates, nltime, nlevel, nlat, 
  #              nlon)
  #   sdatesin: Start date of the input matrix 'YYYYMMDD'.
  #   sdatesout: List of start dates of the output matrix c('YYYYMMDD',
  #                                                         'YYYYMMDD', ...).
  #   nleadtimesout: Number of leadtimes in the output matrix.
  #
  # Returns:
  #   A matrix with the same number of dimensions as the input one, the same 
  #   dimensions 1 and 2 and potentially the same dimensions 5 to 7. Dimensions 
  #   3 and 4 are set by the arguments sdatesout and nleadtimesout.
  #
  # History:
  #   1.0  #  2012-11  (V. Guemas, vguemas@ic3.cat)  #  Original code
  
  #
  #  Input parameters
  # ~~~~~~~~~~~~~~~~~~
  #
  ndims <- length(dim(varin))
  varin <- Enlarge(varin, 7)
  outvar <- array(dim = c(dim(varin)[1:2], length(sdatesout), nleadtimesout,
                  dim(varin)[5:7]))
  yearsin <- as.integer(substr(sdatesin, 1, 4))
  yearsout <- as.integer(substr(sdatesout, 1, 4))
  monthin <- as.integer(substr(sdatesin, 5, 6))
  monthout <- as.integer(substr(sdatesout, 5, 6))
  #
  #  Re-organization
  # ~~~~~~~~~~~~~~~~~
  #
  for (indyear in 1:length(sdatesout)) {
    difmonths <- (yearsin[1] - yearsout[indyear]
                 ) * 12 + monthin[1] - monthout[1]
    if ((difmonths < nleadtimesout) & ((-difmonths) < dim(varin)[4])) {
      outvar[, , indyear, max(difmonths + 1, 1):min(dim(varin)[4] + difmonths,
             nleadtimesout), , , ] <- varin[, , 1, 
                                        max(1 - difmonths, 
                                        1):min(nleadtimesout - difmonths, 
                                        dim(varin)[4]), , , ]
    }
  }
  #
  #  Outputs
  # ~~~~~~~~~
  #
  outvar <- array(outvar, dim = dim(outvar)[1:ndims])
}
