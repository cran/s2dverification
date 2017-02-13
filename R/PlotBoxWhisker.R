PlotBoxWhisker <- function(exp, obs, toptitle = '', ytitle = '', monini = 1, 
                           yearini = 0, freq = 1, expname = "exp 1", 
                           obsname = "obs 1", drawleg = TRUE,
                           fileout = "output_PlotBoxWhisker.ps", 
                           width = 8, height = 5, size_units = 'in', res = 100, ...) {

  # Process the user graphical parameters that may be passed in the call
  ## Graphical parameters to exclude
  excludedArgs <- c("adj", "bty", "cex", "cex.axis", "cex.main", "col", "din", "fin", "lab", "las", "lty", "lwd", "mai", "mar", "mgp", "new", "pch", "ps")
  userArgs <- .FilterUserGraphicArgs(excludedArgs, ...)

  # If there is any filenames to store the graphics, process them
  # to select the right device 
  if (!is.null(fileout)) {
    deviceInfo <- .SelectDevice(fileout = fileout, width = width, height = height, units = size_units, res = res)
    saveToFile <- deviceInfo$fun
    fileout <- deviceInfo$files
  }

  # Checking exp
  if (is.numeric(exp)) {
    if (is.null(dim(exp)) || length(dim(exp)) == 1) {
      dim(exp) <- c(1, length(exp))
    }
  }
  if (!is.numeric(exp) || length(dim(exp)) != 2) {
    stop("Parameter 'exp' must be a numeric vector or array of dimensions c(forecast horizons/start dates) or c(ensemble members, forecast horizons/start dates)")
  }

  # Checking obs
  if (is.numeric(obs)) {
    if (is.null(dim(obs)) || length(dim(obs)) == 1) {
      dim(obs) <- c(1, length(obs))
    }
  }
  if (!is.numeric(obs) || length(dim(obs)) != 2) {
    stop("Parameter 'obs' must be a numeric vector or array of dimensions c(forecast horizons/start dates) or c(1, forecast horizons/start dates)")
  }

  # Checking consistency in exp and obs
  if (dim(exp)[2] != dim(obs)[2]) {
    stop("'exp' and 'obs' must have data for the same amount of time steps.")
  }

  if (!is.character(toptitle) || !is.character(ytitle)) {
    stop("Parameters 'ytitle' and 'toptitle' must be character strings.")
  }

  if (!is.numeric(monini)) {
    stop("'monini' must be a month number, from 1 to 12.")
  }
  if (monini < 1 || monini > 12) {
    stop("'monini' must be >= 1 and <= 12.")
  }

  if (!is.numeric(yearini)) {
    stop("'yearini' must be a month number, from 1 to 12.")
  }

  if (!is.numeric(freq)) {
    stop("'freq' must be a number <= 12.")
  }

  if (!is.character(expname) || !is.character(obsname)) {
    stop("'expname' and 'obsname' must be character strings.")
  }

  if (!is.logical(drawleg)) {
    stop("Parameter 'drawleg' must be either TRUE or FALSE.")
  }

  if (!is.character(fileout) && !is.null(fileout)) {
    stop("Parameter 'fileout' must be a character string.")
  }

  ntimesteps <- dim(exp)[2]
  lastyear <- (monini + (ntimesteps - 1) * 12 / freq - 1) %/% 12 + yearini
  lastmonth <- (monini + (ntimesteps - 1) * 12 / freq - 1) %% 12 + 1
  #
  #  Define some plot parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  labind <- seq(1, ntimesteps)
  months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
              "Oct", "Nov", "Dec")
  labyear <- ((labind - 1) * 12 / freq + monini - 1) %/% 12 + yearini
  labmonth <- months[((labind - 1) * 12 / freq + monini - 1) %% 12 + 1]
  for (jx in 1:length(labmonth)) {
    y2o3dig <- paste("0", as.character(labyear[jx]), sep = "")
    labmonth[jx] <- paste(labmonth[jx], "\nYr ", substr(y2o3dig,
                          nchar(y2o3dig) - 1, nchar(y2o3dig)), sep = "")
  }

  # Open connection to graphical device
  if (!is.null(fileout)) {
    saveToFile(fileout)
  } else if (names(dev.cur()) == 'null device') {
    dev.new(units = size_units, res = res, width = width, height = height)
  }

  # Load the user parameters
  par(userArgs)

  ## Observed time series.
  #pc.o <- ts(obs[1, ], deltat = 1, start = yr1, end = yr2)
  pc.o <- obs[1, ]
  ## Normalization of obs, forecast members. Fabian
  ## Normalization of forecast should be according to ensemble
  ## mean, to keep info on ensemble spread, no? Lauriane pc.o <-
  ## pc.o/sd(pc.o) sd.fc <- apply(exp,c(1),sd)
  ## exp <- exp/sd.fc mn.fc <-
  ## apply(exp,2, mean) exp <-
  ## exp/sd(mn.fc) Produce plot.
  par(mar = c(5, 6, 4, 2))
  boxplot(exp, add = FALSE, main = toptitle, 
    ylab = "", xlab = "", col = "red", lwd = 2, t = "b", 
    axes = FALSE, cex.main = 2, ylim = c(-max(abs(c(exp, pc.o))), max(abs(c(exp, pc.o)))))
  lines(1:ntimesteps, pc.o, lwd = 3, col = "blue")
  abline(h = 0, lty = 1)
  if (drawleg) {
    legend("bottomleft", c(obsname, expname), lty = c(1, 1), lwd = c(3, 
      3), pch = c(NA, NA), col = c("blue", "red"), horiz = FALSE, 
      bty = "n", inset = 0.05)
  }
  ##mtext(1, line = 3, text = tar, cex = 1.9)
  mtext(3, line = -2, text = paste(" AC =", round(cor(pc.o, 
        apply(exp, c(2), mean)), 2)), cex = 1.9, adj = 0)
  axis(2, cex.axis = 2)
  mtext(2, line = 3, text = ytitle, cex = 1.9)
  par(mgp = c(0, 4, 0))
  ##axis(1, c(1:ntimesteps), NA, cex.axis = 2)
  axis(1, seq(1, ntimesteps, by = 1), labmonth, cex.axis = 2)
  box()

  # If the graphic was saved to file, close the connection with the device
  if(!is.null(fileout)) dev.off()
}

