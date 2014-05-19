Plot2VarsVsLTime <- function(var1, var2, toptitle = '', ytitle = '', monini = 1,
                              freq = 12, nticks = NULL, limits = NULL, listexp =
                              c('exp1', 'exp2', 'exp3'), listvars = c('var1',
                             'var2'), biglab = FALSE, hlines = NULL, leg = TRUE,
                             siglev = FALSE, sizetit = 1, fileout = 
                             'output_plot2varsvsltime.eps', show_conf = TRUE) {
  # Plots two input variables havind the same dimensions in a common plot.
  # One plot for all experiments.
  # Input variables should have dimensions (nexp/nmod, nltime).
  #
  # Args:
  #   var1: Matrix of dimensions (nexp/nmod, nltime).
  #   var2: Matrix of dimensions (nexp/nmod, nltime).
  #   toptitle: Main title, optional.
  #   ytitle: Title of Y-axis, optional.
  #   monini: Starting month between 1 and 12. Default = 1.
  #   freq: 1 = yearly, 12 = monthly, 4 = seasonal, ... Default = 12.
  #   nticks: Number of ticks and labels on the x-axis, optional.
  #   limits: c(lower limit, upper limit): limits of the Y-axis, optional.
  #   listexp: List of experiment names, up to three, optional.
  #   listvars: List of names of input variables, optional.
  #   biglab: TRUE/FALSE for presentation/paper plot. Default = FALSE.
  #   hlines: c(a, b, ...) Add horizontal black lines at Y-positions a, b, ...
  #           Default: NULL. 
  #   leg: TRUE/FALSE if legend should be added or not to the plot. Default = TRUE.
  #   siglev: TRUE/FALSE if significance level should replace confidence interval.
  #           Default = FALSE.
  #   sizetit: Multiplicative factor to change title size, optional.
  #   fileout: Name of output ps file.
  #   show_conf: TRUE/FALSE to show/not confidence intervals for input variables.
  #
  # Returns:
  #   This function returns nothing.
  #
  # History:
  #   1.0  #  2013-03  (I. Andreu-Burillo, 
  #                     isabel.andreu-burillo@ic3.cat)  #  Original code
         
  #
  #  Get some arguments
  # ~~~~~~~~~~~~~~~~~~~~
  #
  nvars <- 2

  if (length(dim(var1)) != length(dim(var2))) { 
    print("the two input variables should have the same dimensions")
    stop()
  }
  if (length(dim(var1)) >= 4) { 
    print("dimensions of input variables should be 3")
    stop()
  }
  nleadtime <- dim(var1)[3]
  nexp <- dim(var1)[1]
  var <- array(dim = c(nvars, nexp, 3, nleadtime))
  for (jvar in 1:nvars) {
    varname <- paste("var", as.character(jvar), sep = "")
    var[jvar, , , ] <- get(varname)
    rm(varname)
  }
  if (is.null(limits) == TRUE) {
    ll <- min(var1, na.rm = TRUE)
    ul <- max(var1, na.rm = TRUE)
    if (biglab) {
      ul <- ul + 0.4 * (ul - ll)
    } else {
      ul <- ul + 0.3 * (ul - ll)
    }
  } else {
    ll <- limits[1]
    ul <- limits[2]
  }
  lastyear <- (monini + (nleadtime - 1) * 12 / freq - 1) %/% 12
  lastmonth <- (monini + (nleadtime - 1) * 12 / freq - 1) %% 12 + 1
  empty_ts <- ts(start = c(0000, (monini - 1) %/% (12 / freq) + 1), 
                 end = c(lastyear, (lastmonth - 1) %/% (12 / freq) + 1), 
                 frequency = freq)
  empty <- array(dim = length(empty_ts))
  #
  #  Define some plot parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (is.null(nticks)) {
    if (biglab) {
      nticks <- 5
    } else {
      nticks <- 10
    }
  }
  labind <- seq(1, nleadtime, max(nleadtime %/% nticks, 1))
  months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
              "Oct", "Nov", "Dec")
  labyear <- ((labind - 1) * 12 / freq + monini - 1) %/% 12
  labmonth <- months[((labind - 1) * 12 / freq + monini - 1) %% 12 + 1]
  for (jx in 1:length(labmonth)) {
    y2o3dig <- paste("0", as.character(labyear[jx]), sep = "")
    labmonth[jx] <- paste(labmonth[jx], "\nYr ", substr(y2o3dig, nchar(y2o3dig)
                    - 1, nchar(y2o3dig)), sep = "")
  }
  color <- c("red1", "dodgerblue1", "green1", "orange1", "lightblue1", 
             "deeppink1", "mediumpurple1", "lightgoldenrod1", "olivedrab1",
             "mediumorchid1")
  type <- c(1, 3)
  if (siglev == TRUE) {
    lines <- c("n", "l", "n")
  }
  else{
    lines <- c("l", "l", "l")
  }
  thickness <- array(dim = c(3))
  thickness[1] <- c(1)
  thickness[2] <- c(8)
  thickness[3] <- thickness[1]
  #
  #  Define plot layout
  # ~~~~~~~~~~~~~~~~~~~~
  #
  postscript(fileout, width = 550, height = 300)
  if (biglab) {
    par(mai = c(1.25, 1.4, 0.5, 1), mgp = c(4, 2.5, 0))
    par(cex = 1.3, cex.lab = 2, cex.axis = 1.8)
    cexmain <- 2.2
    legsize <- 1.5
  } else {
    par(mai = c(1, 1.1, 0.5, 0), mgp = c(3, 1.8, 0))
    par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.1)
    cexmain <- 1.5
    legsize <- 1
  }
  plot(empty, ylim = c(ll, ul), xlab = "Time (months)", ylab = ytitle, 
       main = toptitle, cex.main = cexmain * sizetit, axes = FALSE)
  axis(1, at = labind, labels = labmonth)
  axis(2)
  box()
  if (is.null(hlines) != TRUE) { 
    for (jy in 1:length(hlines)) {
      par(new = TRUE) 
      abline(h = hlines[jy])
    }
  }
  #
  #  Loop on experimental data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #                 
  legendnames <- array(dim = nexp * nvars)
  legendthick <- array(dim = nexp * nvars)
  legendsty <- array(dim = nexp * nvars)
  legendcol <- array(dim = nexp * nvars)
  if (show_conf == TRUE) {
    start_line <- 3
    end_line <- 1
  } else {
    start_line <- 2
    end_line <- 2
  }
  for (jint in seq(start_line, end_line, -1)) {
    ind <- 1
    for (jexp in 1:nexp) {
      for (jvar in 1:nvars) {
        par(new = TRUE)
        plot(var[jvar, jexp, jint, ], type = lines[jint], ylim = c(ll, ul), 
             col = color[jexp], lty = type[jvar], lwd = thickness[jint],
             ylab = "", xlab = "", axes = FALSE)
        legendnames[ind] <- paste(listexp[jexp], listvars[jvar])
        legendthick[ind] <- 2
        legendsty[ind] <- type[jvar]
        legendcol[ind] <- color[jexp]
        ind <- ind + 1
      }
    }
  }
  if (leg) {
    legend(1, ul, legendnames, lty = legendsty, lwd = legendthick,
           col = legendcol, cex = legsize)
  }
  dev.off()
}
