AnimateMap <- function(var, lon, lat, toptitle = rep("", 11), sizetit = 1, 
                       units = "", monini = 1, freq = 12, msk95lev = FALSE, 
                       brks = NULL, cols = NULL, filled.continents = FALSE, 
                       lonmin = 0, lonmax = 360, latmin = -90, latmax = 90, 
                       intlon = 20, intlat = 30, drawleg = TRUE, 
                       subsampleg = 1, colNA = "white", equi = TRUE, 
                       fileout = c("output1_animvsltime.gif", 
                                   "output2_animvsltime.gif", 
                                   "output3_animvsltime.gif"), ...) {
  # Process the user graphical parameters that may be passed in the call
  ## Graphical parameters to exclude
  excludedArgs <- c("bg", "col", "fin", "lab", "lend", "new", "pin", "ps")
  userArgs <- .FilterUserGraphicArgs(excludedArgs, ...)

  ## fileout content with extension for consistency between
  ## functions keeping only filename without extension
  ext <- regmatches(fileout, regexpr("\\.[a-zA-Z0-9]*$", fileout))
  if ((length(ext) != 0) && any(ext != ".gif")) {
    .warning("some or all extensions of the filenames provided in 'fileout' are not 'gif'. The extensions are being converted to 'gif'.")
  }
  fileout <- sub("\\.[a-zA-Z0-9]*$", "", fileout)

  #

  # Check var
  if (!is.numeric(var) || !is.array(var)) {
    stop("Parameter 'var' must be a numeric array.")
  }
  if (length(dim(var)) < 3 || length(dim(var)) > 6) {
    stop("Parameter 'var' must be an array with 3 to 6 dimensions.")
  }
  if (length(dim(var)) == 3) {
    var <- InsertDim(var, posdim = 1, lendim = 1)
  }
  if (length(dim(var)) == 4) {
    var <- InsertDim(var, posdim = 2, lendim = 3)
  }
  if (length(dim(var)) == 5) {
    var <- InsertDim(var, posdim = 2, lendim = 1)
  }

  nleadtime <- dim(var)[4]
  nexp <- dim(var)[1]
  nobs <- dim(var)[2]
  nlat <- dim(var)[5]
  nlon <- dim(var)[6]
  if (length(lon) != nlon | length(lat) != nlat) {
    stop("Inconsistent var dimensions / longitudes +  latitudes")
  }
  colorbar <- clim.palette()
  if (is.null(brks) == TRUE) {
    ll <- signif(min(var[, , 2, , , ], na.rm = TRUE), 4)
    ul <- signif(max(var[, , 2, , , ], na.rm = TRUE), 4)
    if (is.null(cols) == TRUE) {
      cols <- colorbar(10)
    }
    nlev <- length(cols)
    brks <- signif(seq(ll, ul, (ul - ll)/nlev), 4)
  } else {
    if (is.null(cols) == TRUE) {
      nlev <- length(brks) - 1
      cols <- colorbar(nlev)
    } else {
      nlev <- length(cols)
    }
  }
  lon[which(lon < lonmin)] <- lon[which(lon < lonmin)] + 360
  lon[which(lon > lonmax)] <- lon[which(lon > lonmax)] - 360
  latb <- sort(lat[which(lat >= latmin & lat <= latmax)], index.return = TRUE)
  lonb <- sort(lon[which(lon >= lonmin & lon <= lonmax)], index.return = TRUE)
  
  # Define some plot parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  labind <- 1:nleadtime
  months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", 
    "Aug", "Sep", "Oct", "Nov", "Dec")
  years <- ((labind - 1) * 12/freq + monini - 1)%/%12
  suffixtit <- months[((labind - 1) * 12/freq + monini - 1)%%12 + 
    1]
  for (jx in 1:nleadtime) {
    y2o3dig <- paste("0", as.character(years[jx]), sep = "")
    suffixtit[jx] <- paste(suffixtit[jx], "-", substr(y2o3dig, 
      nchar(y2o3dig) - 1, nchar(y2o3dig)), sep = "")
  }
  
  # Loop on experimental & observational data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (jexp in 1:nexp) {
    for (jobs in 1:nobs) {
      postscript(paste(fileout[(jexp - 1) * nobs + jobs], 
        ".png", sep = ""), width = 550, height = 300, 
        bg = "white")
      # Load the user parameters
      par(userArgs)
      for (jt in 1:nleadtime) {
        title <- paste(toptitle[(jexp - 1) * nobs + jobs], 
          " Time=", suffixtit[jt], sep = "")
        varbis <- var[jexp, jobs, 2, jt, which(lat >= 
          latmin & lat <= latmax), which(lon >= lonmin & 
          lon <= lonmax)]
        varbis <- varbis[latb$ix, lonb$ix]
        flag <- array(FALSE, dim(varbis))
        if (msk95lev) {
          flag[which(var[jexp, jobs, 1, jt, latb$ix, 
          lonb$ix] > 0 & var[jexp, jobs, 3, jt, latb$ix, 
          lonb$ix] > 0)] <- TRUE
          flag[which(var[jexp, jobs, 1, jt, latb$ix, 
          lonb$ix] < 0 & var[jexp, jobs, 3, jt, latb$ix, 
          lonb$ix] < 0)] <- TRUE
        }
        varbis[which(varbis <= min(brks))] <- min(brks) + 
          (max(brks) - min(brks))/1000
        varbis[which(varbis >= max(brks))] <- max(brks) - 
          (max(brks) - min(brks))/1000
        if (equi) {
          PlotEquiMap(t(varbis), lonb$x, latb$x, toptitle = title, 
          sizetit = sizetit, units = units, filled.continents = filled.continents, 
          dots = t(flag), brks = brks, cols = cols, 
          intxlon = intlon, intylat = intlat, drawleg = drawleg, 
          subsampleg = subsampleg, colNA = colNA, ...)
        } else {
          PlotStereoMap(t(varbis), lonb$x, latb$x, latlims = c(latmin, 
          latmax), toptitle = title, sizetit = sizetit, 
          units = units, filled.continents = filled.continents, 
          dots = t(flag), brks = brks, cols = cols, 
          intlat = intlat, drawleg = drawleg, subsampleg = subsampleg, 
          colNA = colNA, ...)
        }
      }
      dev.off()
      system(paste("convert -rotate 90 -loop 10 -delay 50 ", 
        fileout[(jexp - 1) * nobs + jobs], ".png ", fileout[(jexp - 
          1) * nobs + jobs], ".gif", sep = ""))
      file.remove(paste0(fileout[(jexp - 1) * nobs + jobs], 
        ".png"))
    }
  }
} 
