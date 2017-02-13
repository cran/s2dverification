Histo2Hindcast <- function(varin, sdatesin, sdatesout, nleadtimesout) {
  #
  #  Input parameters
  # ~~~~~~~~~~~~~~~~~~
  #
  ndims <- length(dim(varin))
  varin <- Enlarge(varin, 7)
  outvar <- array(dim = c(dim(varin)[1:2], length(sdatesout), nleadtimesout,
                  dim(varin)[5:7]))
  yearsin <- as.numeric(substr(sdatesin, 1, 4))
  yearsout <- as.numeric(substr(sdatesout, 1, 4))
  monthin <- as.numeric(substr(sdatesin, 5, 6))
  monthout <- as.numeric(substr(sdatesout, 5, 6))
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
