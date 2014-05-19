ACC <- function(var_exp, var_obs, lon = NULL, lat = NULL, lonlatbox = NULL) {
  # Matrix var_exp & var_obs should have dimensions (nexp/nobs, nsdates, 
  # nltimes, nlat, nlon).
  # ACC computes the Anomaly Correlation Coefficient for each jexp in 
  # 1:nexp and each jobs in 1:nobs which gives nexp x nobs ACC for each 
  # startdate and each leadtime.
  # The confidence interval is computed by a Fisher transformation.
  # The significance level relies on a one-sided student-T distribution.
  # A domain can be selected by providing the list of longitudes/latitudes
  # (lon/lat) of the grid together with the corner of the domain:
  #   lonlatbox = c(lonmin, lonmax, latmin, latmax)
  # 
  # Args:
  #   var_exp: Matrix of experimental data.
  #   var_obs: Matrix of observational data, same dimensions as var_exp except 
  #            along the first dimension.
  #   lon: Array of longitudes of the var_exp/var_obs grids, optional.
  #   lat: Array of latitudes of the var_exp/var_obs grids, optional.
  #   lonlatbox: Domain to select c(lonmin, lonmax, latmin, latmax), optional.
  # 
  # Returns:
  #   $ACC: Matrix with c(nexp, nobs, nsdates, nleadtimes, 4) dimensions.
  #        The fourth dimension of length 4 corresponds to the lower limit of
  #        the 95% confidence interval, the computed ACC, the upper limit of 
  #        the 95% confidence interval and the 95% significance level given by 
  #        a one-sided T-test.
  #   $MACC: Mean Anomaly Correlation Coefficient with c(nexp, nobs, nleadtimes) 
  #         dimensions.
  # 
  # History:
  #   1.0  #  2013-08  (V. Guemas, vguemas@ic3.cat)  #  Original code
  
  dimsvar <- dim(var_exp)
  if (length(dimsvar) != 5) {
    stop("var_exp & var_obs should have dimensions (nexp/nsobs, nsdates, nltimes, nlat, nlon) ")
  }
  for (iind in 2:length(dimsvar)) {
    if (dim(var_obs)[iind] != dimsvar[iind]) {
      stop("var_exp & var_obs must have same dimensions except the first one (number of experiments or number of observational datasets) ")
    }
  }
  nexp <- dimsvar[1]
  nobs <- dim(var_obs)[1]
  nsdates <- dimsvar[2]
  nltimes <- dimsvar[3]
  nlat <- dimsvar[4]
  nlon <- dimsvar[5]
  if (is.null(lon) == FALSE & is.null(lat) == FALSE & 
      is.null(lonlatbox) == FALSE) {
    for (jind in 1:2) {
      while (lonlatbox[jind] < 0) {
        lonlatbox[jind] <- lonlatbox[jind] + 360
      }
      while (lonlatbox[jind] > 360) {
        lonlatbox[jind] <- lonlatbox[jind] - 360
      }
    }
    indlon <- which((lon >= lonlatbox[1] & lon <= lonlatbox[2]) | 
                    (lonlatbox[1] > lonlatbox[2] & (lon > lonlatbox[1] | 
                    lon < lonlatbox[2])))
    indlat <- which(lat >= lonlatbox[3] & lat <= lonlatbox[4])
  } else {
    indlon <- 1:nlon
    indlat <- 1:nlat
  }
  ACC <- array(dim = c(nexp, nobs, nsdates, nltimes, 4))
  MACC <- array(0, dim = c(nexp, nobs, nltimes))
  for (jexp in 1:nexp) {
    for (jobs in 1:nobs) {
      for (jltime in 1:nltimes) {
        top <- 0
        bottom1 <- 0
        bottom2 <- 0
        numdates <- 0
        for (jdate in 1:nsdates) {
          tmp1 <- array(var_exp[jexp, jdate, jltime, indlat, indlon], 
                        dim = length(indlon) * length(indlat))
          tmp2 <- array(var_obs[jobs, jdate, jltime, indlat, indlon], 
                        dim = length(indlon) * length(indlat))
          notna <- which(is.na(tmp1) == FALSE & is.na(tmp2) == FALSE)
          if (length(sort(tmp1)) > 0 & length(sort(tmp2)) > 2) {
            toto <- sum(tmp1[notna] * tmp2[notna]) / sqrt((sum(
                    tmp1[notna] ^ 2)) * (sum(tmp2[notna] ^ 2)))
            ACC[jexp, jobs, jdate, jltime, 2] <- toto
            eno <- Eno(tmp2, 1)
            t <- qt(0.95, eno - 2)
            ACC[jexp, jobs, jdate, jltime, 4] <- sqrt((t * t) / (
                                                 (t * t) + eno - 2))
            ACC[jexp, jobs, jdate, jltime, 1] <- tanh(atanh(toto) + qnorm(
                                                 0.975) / sqrt(eno - 3))
            ACC[jexp, jobs, jdate, jltime, 3] <- tanh(atanh(toto) + qnorm(
                                                 0.025) / sqrt(eno - 3))
            top <- top + sum(tmp1[notna] * tmp2[notna])
            bottom1 <- bottom1 + sum(tmp1[notna] ^ 2)
            bottom2 <- bottom2 + sum(tmp2[notna] ^ 2)
            numdates <- numdates + 1
          }
        }
        if (numdates > 0) {
          MACC[jexp, jobs, jltime] <- top / sqrt(bottom1 * bottom2)
        }
      }
    }
  }
  invisible(list(ACC = ACC, MACC = MACC))
}
