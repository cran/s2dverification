PlotEquiMap <- function(var, lon, lat, varu = NULL, varv = NULL,
                        toptitle = NULL, sizetit = NULL, units = NULL, 
                        brks = NULL, cols = NULL, bar_limits = NULL, 
                        triangle_ends = NULL, col_inf = NULL, col_sup = NULL, 
                        colNA = NULL, color_fun = clim.palette(),
                        square = TRUE, filled.continents = NULL,
                        coast_color = NULL, coast_width = 1,
                        contours = NULL, brks2 = NULL, contour_lwd = 0.5,
                        contour_color = 'black', contour_lty = 1,
                        contour_label_scale = 1,
                        dots = NULL, dot_symbol = 4, dot_size = 1, 
                        arr_subsamp = floor(length(lon) / 30), arr_scale = 1, 
                        arr_ref_len = 15, arr_units = "m/s", 
                        arr_scale_shaft = 1, arr_scale_shaft_angle = 1,
                        axelab = TRUE, labW = FALSE, 
                        intylat = 20, intxlon = 20, 
                        axes_tick_scale = 1, axes_label_scale = 1,
                        drawleg = TRUE, subsampleg = NULL, 
                        bar_extra_labels = NULL, draw_bar_ticks = TRUE, 
                        draw_separators = FALSE, triangle_ends_scale = 1, 
                        bar_label_digits = 4, bar_label_scale = 1, 
                        units_scale = 1, bar_tick_scale = 1, 
                        bar_extra_margin = rep(0, 4),
                        boxlim = NULL, boxcol = 'purple2', boxlwd = 5, 
                        margin_scale = rep(1, 4), title_scale = 1, 
                        numbfig = NULL, fileout = NULL, 
                        width = 8, height = 5, size_units = 'in', 
                        res = 100, ...) {
  # Process the user graphical parameters that may be passed in the call
  ## Graphical parameters to exclude
  excludedArgs <- c("cex", "cex.axis", "cex.lab", "cex.main", "col", "din", "fig", "fin", "lab", "las", "lty", "lwd", "mai", "mar", "mgp", "new", "oma", "ps", "tck")
  userArgs <- .FilterUserGraphicArgs(excludedArgs, ...)

  # If there is any filenames to store the graphics, process them
  # to select the right device 
  if (!is.null(fileout)) {
    deviceInfo <- .SelectDevice(fileout = fileout, width = width, height = height, units = size_units, res = res)
    saveToFile <- deviceInfo$fun
    fileout <- deviceInfo$files
  }

  # Preliminar check of dots, contours, varu, varv, lon, lat
  if (!is.null(dots)) {
    if (!is.array(dots) || !(length(dim(dots)) %in% c(2, 3))) {
      stop("Parameter 'dots' must be a logical array with two or three dimensions.")
    }
    if (length(dim(dots)) == 2) {
      dim(dots) <- c(1, dim(dots))
    }
  }
  if (!is.null(contours)) {
    if (!is.array(contours) || !(length(dim(contours)) == 2)) {
      stop("Parameter 'contours' must be a numerical array with two dimensions.")
    }
  }
  if (!is.null(varu) && !is.null(varv)) {
    if (!is.array(varu) || !(length(dim(varu)) == 2)) {
      stop("Parameter 'varu' must be a numerical array with two dimensions.")
    }
    if (!is.array(varv) || !(length(dim(varv)) == 2)) {
      stop("Parameter 'varv' must be a numerical array with two dimensions.")
    }
  } else if (!is.null(varu) || !is.null(varv)) {
    stop("Only one of the components 'varu' or 'varv' has been provided. Both must be provided.")
  }
  if (!is.numeric(lon) || !is.numeric(lat)) {
    stop("Parameters 'lon' and 'lat' must be numeric vectors.")
  }

  # Check var
  if (!is.array(var)) {
    stop("Parameter 'var' must be a numeric array.")
  }
  if (length(dim(var)) > 2) {
    var <- drop(var)
    dim(var) <- head(c(dim(var), 1, 1), 2)
  }
  if (length(dim(var)) > 2) {
    stop("Parameter 'var' must be a numeric array with two dimensions. See PlotMultiMap() for multi-pannel maps or AnimateMap() for animated maps.")
  } else if (length(dim(var)) < 2) {
    stop("Parameter 'var' must be a numeric array with two dimensions.")
  }
  dims <- dim(var)
  # Transpose the input matrices because the base plot functions work directly 
  # with dimensions c(lon, lat).
  if (dims[1] != length(lon) || dims[2] != length(lat)) {
    if (dims[1] == length(lat) && dims[2] == length(lon)) {
      var <- t(var)
      if (!is.null(varu)) varu <- t(varu)
      if (!is.null(varv)) varv <- t(varv)
      if (!is.null(contours)) contours <- t(contours)
      if (!is.null(dots)) dots <- aperm(dots, c(1, 3, 2))
      dims <- dim(var)
    }
  }

  # Check lon
  if (length(lon) != dims[1]) {
    stop("Parameter 'lon' must have as many elements as the number of cells along longitudes in the input array 'var'.")
  }

  # Check lat
  if (length(lat) != dims[2]) {
    stop("Parameter 'lat' must have as many elements as the number of cells along longitudes in the input array 'var'.")
  }

  # Check varu and varv
  if (!is.null(varu) && !is.null(varv)) {
    if (dim(varu)[1] != dims[1] || dim(varu)[2] != dims[2]) {
      stop("Parameter 'varu' must have same number of longitudes and latitudes as 'var'.")
    }
    if (dim(varv)[1] != dims[1] || dim(varv)[2] != dims[2]) {
      stop("Parameter 'varv' must have same number of longitudes and latitudes as 'var'.")
    }
  }

  # Check toptitle
  if (is.null(toptitle) || is.na(toptitle)) {
    toptitle <- ''
  }
  if (!is.character(toptitle)) {
    stop("Parameter 'toptitle' must be a character string.")
  }

  # Check sizetit
  if (!is.null(sizetit)) {
    .warning("Parameter 'sizetit' is obsolete. Use 'title_scale' instead.")
    if (!is.numeric(sizetit) || length(sizetit) != 1) {
      stop("Parameter 'sizetit' must be a single numeric value.")
    }
    title_scale <- sizetit
  }

  var_limits <- c(min(var, na.rm = TRUE), max(var, na.rm = TRUE))
  # Check: brks, cols, subsampleg, bar_limits, color_fun, bar_extra_labels, draw_bar_ticks
  #        draw_separators, triangle_ends_scale, label_scale, units, units_scale, 
  #        bar_label_digits
  # Build: brks, cols, bar_limits, col_inf, col_sup
  colorbar <- ColorBar(brks, cols, FALSE, subsampleg, bar_limits, var_limits,
                       triangle_ends, col_inf, col_sup, color_fun, FALSE, 
                       extra_labels = bar_extra_labels, draw_ticks = draw_bar_ticks,
                       draw_separators = draw_separators, 
                       triangle_ends_scale = triangle_ends_scale,
                       label_scale = bar_label_scale, title = units, 
                       title_scale = units_scale, tick_scale = bar_tick_scale,
                       extra_margin = bar_extra_margin, label_digits = bar_label_digits)
  brks <- colorbar$brks
  cols <- colorbar$cols
  col_inf <- colorbar$col_inf
  col_sup <- colorbar$col_sup
  bar_limits <- c(head(brks, 1), tail(brks, 1))

  # Check colNA
  if (is.null(colNA)) {
    if ('na_color' %in% names(attributes(cols))) {
      colNA <- attr(cols, 'na_color')
      if (!.IsColor(colNA)) {
        stop("The 'na_color' provided as attribute of the colour vector must be a valid colour identifier.")
      }
    } else {
      colNA <- 'pink'
    }
  } else if (!.IsColor(colNA)) {
    stop("Parameter 'colNA' must be a valid colour identifier.")
  }

  # Check square
  if (!is.logical(square)) {
    stop("Parameter 'square' must be logical.")
  }

  # Check filled.continents
  if (is.null(filled.continents)) {
    if (!square) {
      filled.continents <- FALSE
    } else {
      filled.continents <- TRUE
    }
  }
  if (!.IsColor(filled.continents) && !is.logical(filled.continents)) {
    stop("Parameter 'filled.continents' must be logical or a colour identifier.")
  } else if (!is.logical(filled.continents)) {
    continent_color <- filled.continents
    filled.continents <- TRUE
  } else if (filled.continents) {
    continent_color <- gray(0.5)
  }

  # Check coast_color
  if (is.null(coast_color)) {
    if (filled.continents) {
      coast_color <- continent_color
    } else {
      coast_color <- 'black'
    }
  }
  if (!.IsColor(coast_color)) {
    stop("Parameter 'coast_color' must be a valid colour identifier.")
  }

  # Check coast_width
  if (!is.numeric(coast_width)) {
    stop("Parameter 'coast_width' must be numeric.")
  }

  # Check contours
  if (!is.null(contours)) {
    if (dim(contours)[1] != dims[1] || dim(contours)[2] != dims[2]) {
      stop("Parameter 'contours' must have the same number of longitudes and latitudes as 'var'.")
    }
  }

  # Check brks2
  if (is.null(brks2)) {
    if (is.null(contours)) { 
      if (!square) {
        brks2 <- brks
        contours <- var 
      }
    } else {
      ll <- signif(min(contours, na.rm = TRUE), 2)
      ul <- signif(max(contours, na.rm = TRUE), 2)
      brks2 <- signif(seq(ll, ul, length.out = length(brks)), 2)
    }
  }

  # Check contour_lwd
  if (!is.numeric(contour_lwd)) {
    stop("Parameter 'contour_lwd' must be numeric.")
  }

  # Check contour_color
  if (!.IsColor(contour_color)) {
    stop("Parameter 'contour_color' must be a valid colour identifier.")
  }

  # Check contour_lty
  if (!is.numeric(contour_lty) && !is.character(contour_lty)) {
    stop("Parameter 'contour_lty' must be either a number or a character string.")
  }

  # Check contour_label_scale
  if (!is.numeric(contour_label_scale)) {
    stop("Parameter 'contour_label_scale' must be numeric.")
  }

  # Check dots, dot_symbol and dot_size
  if (!is.null(dots)) {
    if (dim(dots)[2] != dims[1] || dim(dots)[3] != dims[2]) {
      stop("Parameter 'dots' must have the same number of longitudes and latitudes as 'var'.")
    }
    if (!is.numeric(dot_symbol) && !is.character(dot_symbol)) {
      stop("Parameter 'dot_symbol' must be a numeric or character string vector.")
    }
    if (length(dot_symbol) == 1) {
      dot_symbol <- rep(dot_symbol, dim(dots)[1])
    } else if (length(dot_symbol) < dim(dots)[1]) {
      stop("Parameter 'dot_symbol' does not contain enough symbols.")
    }
    if (!is.numeric(dot_size)) {
      stop("Parameter 'dot_size' must be numeric.")
    }
    if (length(dot_size) == 1) {
      dot_size <- rep(dot_size, dim(dots)[1])
    } else if (length(dot_size) < dim(dots)[1]) {
      stop("Parameter 'dot_size' does not contain enough sizes.")
    }
  }

  # Check arrow parameters
  if (!is.numeric(arr_subsamp)) {
    stop("Parameter 'arr_subsamp' must be numeric.")
  }
  if (!is.numeric(arr_scale)) {
    stop("Parameter 'arr_scale' must be numeric.")
  }
  if (!is.numeric(arr_ref_len)) {
    stop("Parameter 'arr_ref_len' must be numeric.")
  }
  if (!is.character(arr_units)) {
    stop("Parameter 'arr_units' must be character.")
  }
  if (!is.numeric(arr_scale_shaft)) {
    stop("Parameter 'arr_scale_shaft' must be numeric.")
  }
  if (!is.numeric(arr_scale_shaft_angle)) {
    stop("Parameter 'arr_scale_shaft_angle' must be numeric.")
  }

  # Check axis parameters
  if (!is.logical(axelab)) {
    stop("Parameter 'axelab' must be logical.")
  }
  if (!is.logical(labW)) {
    stop("Parameter 'labW' must be logical.")
  }
  if (!is.numeric(intylat)) {
    stop("Parameter 'intylat' must be numeric.")
  } else {
    intylat <- round(intylat)
  }
  if (!is.numeric(intxlon)) {
    stop("Parameter 'intxlon' must be numeric.")
  } else {
    intxlon <- round(intxlon)
  }

  # Check legend parameters
  if (!is.logical(drawleg)) {
    stop("Parameter 'drawleg' must be logical.")
  }

  # Check box parameters
  if (!is.null(boxlim)) {
    if (!is.list(boxlim)) {
      boxlim <- list(boxlim)
    }
    for (i in 1:length(boxlim)) {
      if (!is.numeric(boxlim[[i]]) || length(boxlim[[i]]) != 4) {
        stop("Parameter 'boxlim' must be a a numeric vector or a list of numeric vectors of length 4 (with W, S, E, N box limits).")
      }
    }
    if (!is.character(boxcol)) {
      stop("Parameter 'boxcol' must be a character string or a vector of character strings.")
    } else {
      if (length(boxlim) != length(boxcol)) {
        if (length(boxcol) == 1) {
          boxcol <- rep(boxcol, length(boxlim))
        } else {
          stop("Parameter 'boxcol' must have a colour for each box in 'boxlim' or a single colour for all boxes.")
        }
      }
    }
    if (!is.numeric(boxlwd)) {
      stop("Parameter 'boxlwd' must be numeric.")
    } else {
      if (length(boxlim) != length(boxlwd)) {
        if (length(boxlwd) == 1) {
          boxlwd <- rep(boxlwd, length(boxlim))
        } else {
          stop("Parameter 'boxlwd' must have a line width for each box in 'boxlim' or a single line width for all boxes.")
        }
      }
    }
  }

  # Check margin_scale
  if (!is.numeric(margin_scale) || length(margin_scale) != 4) {
    stop("Parameter 'margin_scale' must be a numeric vector of length 4.")
  }

  # Check title_scale
  if (!is.numeric(title_scale)) {
    stop("Parameter 'title_scale' must be numeric.")
  }

  # Check axes_tick_scale
  if (!is.numeric(axes_tick_scale)) {
    stop("Parameter 'axes_tick_scale' must be numeric.")
  }

  # Check axes_label_scale
  if (!is.numeric(axes_label_scale)) {
    stop("Parameter 'axes_label_scale' must be numeric.")
  }

  # Check numbfig
  if (!is.null(numbfig)) {
    if (!is.numeric(numbfig)) {
      stop("Parameter 'numbfig' must be numeric.")
    } else {
      numbfig <- round(numbfig)
      scale <- 1 / numbfig ** 0.3
      axes_tick_scale <- axes_tick_scale * scale
      axes_label_scale <- axes_label_scale * scale
      title_scale <- title_scale * scale
      margin_scale <- margin_scale * scale
      arr_scale <- arr_scale * scale
      dot_size <- dot_size * scale
      contour_label_scale <- contour_label_scale * scale
      contour_lwd <- contour_lwd * scale
    }
  }

  #library(GEOmap)
  #library(geomapdata)
  #library(maps)
  utils::data(coastmap, package = 'GEOmap', envir = environment())
  #
  #  Input arguments 
  # ~~~~~~~~~~~~~~~~~
  #
  latb <- sort(lat, index.return = TRUE)
  dlon <- lon[2:dims[1]] - lon[1:(dims[1] - 1)]
  wher <- which(dlon > (mean(dlon) + 1))
  if (length(wher) > 0) {
    lon[(wher + 1):dims[1]] <- lon[(wher + 1):dims[1]] - 360
  }
  lonb <- sort(lon, index.return = TRUE)
  latmin <- floor(min(lat) / 10) * 10
  latmax <- ceiling(max(lat) / 10) * 10
  lonmin <- floor(min(lon) / 10) * 10
  lonmax <- ceiling(max(lon) / 10) * 10
  if (min(lon) < 0) {
    continents <- 'world'
  } else {
    continents <- 'world2'
  }

  #
  #  Plotting the map
  # ~~~~~~~~~~~~~~~~~~
  #

  # Open connection to graphical device
  if (!is.null(fileout)) {
    saveToFile(fileout)
  } else if (names(dev.cur()) == 'null device') {
    dev.new(units = size_units, res = res, width = width, height = height)
  }

  #
  #  Defining the layout
  # ~~~~~~~~~~~~~~~~~~~~~
  # 
  if (drawleg) {
    margin_scale[1] <- margin_scale[1] - 1
  }
  margins <- rep(0.4, 4) * margin_scale
  cex_title <- 2 * title_scale
  cex_axes_labels <- 1.3 * axes_label_scale
  cex_axes_ticks <- -0.5 * axes_tick_scale
  spaceticklab <- 0
  if (axelab) { 
    ypos <- seq(latmin, latmax, intylat)
    xpos <- seq(lonmin, lonmax, intxlon)
    letters <- array('', length(ypos))
    letters[ypos < 0] <- 'S'
    letters[ypos > 0] <- 'N'
    ylabs <- paste(as.character(abs(ypos)), letters, sep = '')
    letters <- array('', length(xpos))
    if (labW) {
      xpos2 <- xpos  
      xpos2[xpos2 > 180] <- 360 - xpos2[xpos2 > 180]
    }
    letters[xpos < 0] <- 'W'
    letters[xpos > 0] <- 'E'
    if (labW) {
      letters[xpos == 0] <- ' '
      letters[xpos == 180] <- ' '
      letters[xpos > 180] <- 'W'  
      xlabs <- paste(as.character(abs(xpos2)), letters, sep = '')
    } else {
      xlabs <- paste(as.character(abs(xpos)), letters, sep = '')
    }
    spaceticklab <- max(-cex_axes_ticks, 0)
    margins[1] <- margins[1] + 1.2 * cex_axes_labels + spaceticklab
    margins[2] <- margins[2] + 1.2 * cex_axes_labels + spaceticklab
  }
  bar_extra_margin[2] <- bar_extra_margin[2] + margins[2]
  bar_extra_margin[4] <- bar_extra_margin[4] + margins[4]
  if (toptitle != '') {
    margins[3] <- margins[3] + cex_title + 1
  }
  if (!is.null(varu)) {
    margins[1] <- margins[1] + 2.2 * units_scale
  }

  if (drawleg) {
    layout(matrix(1:2, ncol = 1, nrow = 2), heights = c(5, 1))
  }
  plot.new()
  # Load the user parameters
  par(userArgs)
  par(mar = margins, cex.main = cex_title, cex.axis = cex_axes_labels,
      mgp = c(0, spaceticklab, 0), las = 0)
  plot.window(xlim = range(lonb$x, finite = TRUE), ylim = range(latb$x, finite = TRUE),
              xaxs = 'i', yaxs = 'i')
  if (axelab) {
    axis(2, at = ypos, labels = ylabs, cex.axis = cex_axes_labels, tcl = cex_axes_ticks,
         mgp = c(0, spaceticklab + 0.2, 0))
    axis(1, at = xpos, labels = xlabs, cex.axis = cex_axes_labels, tcl = cex_axes_ticks,
         mgp = c(0, spaceticklab + cex_axes_labels / 2 - 0.3, 0))
  }
  title(toptitle, cex.main = cex_title)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colNA)
  col_inf_image <- ifelse(is.null(col_inf), colNA, col_inf)
  col_sup_image <- ifelse(is.null(col_sup), colNA, col_sup)
  if (square) {
    image(lonb$x, latb$x, var[lonb$ix, latb$ix], col = c(col_inf_image, cols, col_sup_image), 
          breaks = c(-.Machine$double.xmax, brks, .Machine$double.xmax),
          axes = FALSE, xlab = "", ylab = "", add = TRUE)          
  } else {
    .filled.contour(lonb$x, latb$x, var[lonb$ix, latb$ix], 
                    levels = c(.Machine$double.xmin, brks, .Machine$double.xmax), 
                    col = c(col_inf_image, cols, col_sup_image))
  }
  if (!is.null(contours)) {
    contour(lonb$x, latb$x, contours[lonb$ix, latb$ix], levels = brks2,
            method = "edge", add = TRUE, 
            labcex = cex_axes_labels, lwd = contour_lwd, lty = contour_lty,
            col = contour_color)
  }

  #
  #  Adding black dots or symbols
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (!is.null(dots)) {
    data_avail <- !is.na(var)
    for (counter in 1:(dim(dots)[1])) {
      points <- which(dots[counter, , ] & data_avail, arr.ind = TRUE)
      points(lon[points[, 1]], lat[points[, 2]], 
             pch = dot_symbol[counter], 
             cex = dot_size[counter] * 3 / sqrt(sqrt(length(var))),
             lwd = dot_size[counter] * 3 / sqrt(sqrt(length(var))))
    }
  }
  #
  #  Plotting continents
  # ~~~~~~~~~~~~~~~~~~~~~
  # 
  coast <- map(continents, interior = FALSE, wrap = TRUE,
               fill = filled.continents, add = TRUE, plot = FALSE)
  if (filled.continents) {
    old_lwd <- par('lwd')
    par(lwd = coast_width)
    if (min(lon) >= 0) {
      ylat <- latmin:latmax
      xlon <- lonmin:lonmax
      proj <- GEOmap::setPROJ(1, LON0 = mean(xlon),
                              LAT0 = mean(ylat), LATS = ylat, LONS = xlon)
      lakes <- which(coastmap$STROKES$col == "blue")
      coastmap$STROKES$col[which(coastmap$STROKES$col != "blue")] <- continent_color
      coastmap$STROKES$col[lakes] <- "white"
      par(new = TRUE)
      GEOmap::plotGEOmap(coastmap, PROJ = proj, border = coast_color, 
                         add = TRUE, lwd = coast_width)
    } else {
      polygon(coast, col = continent_color, border = coast_color, lwd = coast_width)
    }
    par(lwd = old_lwd)
  } else {
    lines(coast, col = coast_color, lwd = coast_width)
  }
  box()
  # Draw rectangle on the map
  if (!is.null(boxlim)) {
    counter <- 1
    for (box in boxlim) {
      if (box[1] > box[3]) {
        box[1] <- box[1] - 360
      }
      if (length(box) != 4) {
        stop(paste("The", counter, "st box defined in the parameter 'boxlim' is ill defined."))
      } else if (box[2] < latmin || box[4] > latmax || 
                 box[1] < lonmin || box[3] > lonmax) {
        stop(paste("The limits of the", counter, "st box defined in the parameter 'boxlim' are invalid."))
      } else if (box[1] < 0 && box[3] > 0) {
        #segments south
        segments(box[1], box[2], 0, box[2], col = boxcol[counter], lwd = boxlwd[counter])
        segments(0, box[2], box[3], box[2], col = boxcol[counter], lwd = boxlwd[counter]) 
        #segments north
        segments(box[1], box[4], 0, box[4], col = boxcol[counter], lwd = boxlwd[counter])
        segments(0, box[4], box[3], box[4], col = boxcol[counter], lwd = boxlwd[counter]) 
        #segments west
        segments(box[1], box[2], box[1], box[4], col = boxcol[counter], 
                 lwd = boxlwd[counter])  
        #segments est
        segments(box[3], box[2], box[3],box[4], col = boxcol[counter], 
                 lwd = boxlwd[counter])          
      } else {
        rect(box[1], box[2], box[3], box[4], border = boxcol[counter], col = NULL, 
             lwd = boxlwd[counter], lty = 'solid')
      }
      counter <- counter + 1
    }
  }
  #
  #  PlotWind
  # ~~~~~~~~~~
  #
  if (!is.null(varu) && !is.null(varv)) {
    # Create a two dimention array of longitude and latitude
    lontab <- InsertDim(lonb$x, 2, length( latb$x))
    lattab <- InsertDim(latb$x, 1, length( lonb$x))
    varplotu <- varu[lonb$ix, latb$ix]
    varplotv <- varv[lonb$ix, latb$ix]

    # Select a subsample af the points to an arrow
    #for each "subsample" grid point
    sublon <- seq(1,length(lon), arr_subsamp)
    sublat <- seq(1,length(lat), arr_subsamp)

    uaux <- lontab[sublon, sublat] + varplotu[sublon, sublat] * 0.5 * arr_scale
    vaux <- lattab[sublon, sublat] + varplotv[sublon, sublat] * 0.5 * arr_scale

    lenshaft <- 0.18 * arr_scale * arr_scale_shaft
    angleshaft <- 12 * arr_scale_shaft_angle
    # Plot Wind
    arrows(lontab[sublon, sublat], lattab[sublon, sublat],
           uaux, vaux,
           angle = angleshaft,
           length = lenshaft)
    
    # Plotting an arrow at the bottom of the plot for the legend
    posarlon <- lonb$x[1] + (lonmax - lonmin) * 0.1
    posarlat <- latmin - ((latmax - latmin) + 1) / par('pin')[2] * 
                         (spaceticklab + 0.2 + cex_axes_labels + 0.6 * units_scale) * par('csi')

    arrows(posarlon, posarlat,
           posarlon + 0.5 * arr_scale * arr_ref_len, posarlat,
           length = lenshaft, angle = angleshaft,
           xpd = TRUE)
    #save the parameter value
    xpdsave <- par('xpd')
    #desactivate xpd to be able to plot in margen
    par(xpd = NA)
    #plot text
    mtext(paste(as.character(arr_ref_len), arr_units, sep = ""),
          line = spaceticklab + 0.2 + cex_axes_labels + 1.2 * units_scale, side = 1,
          at = posarlon + (0.5 * arr_scale * arr_ref_len) / 2,
          cex = units_scale)
    #come back to the previous xpd value
    par(xpd = xpdsave)
  }
  #
  #  Colorbar
  # ~~~~~~~~~~
  #
  if (drawleg) {
    ColorBar(brks, cols, FALSE, subsampleg, bar_limits, var_limits, 
             triangle_ends, col_inf = col_inf, col_sup = col_sup, 
             extra_labels = bar_extra_labels, draw_ticks = draw_bar_ticks,
             draw_separators = draw_separators, title = units, 
             title_scale = units_scale, triangle_ends_scale = triangle_ends_scale, 
             label_scale = bar_label_scale, tick_scale = bar_tick_scale,
             extra_margin = bar_extra_margin, label_digits = bar_label_digits)
  }

  # If the graphic was saved to file, close the connection with the device
  if (!is.null(fileout)) dev.off()

  invisible(list(brks = brks, cols = cols, col_inf = col_inf, col_sup = col_sup))
}
