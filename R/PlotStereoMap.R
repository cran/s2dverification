PlotStereoMap <- function(var, lon, lat, latlims = c(60, 90), 
                          toptitle = NULL, sizetit = NULL, units = NULL, 
                          brks = NULL, cols = NULL, bar_limits = NULL, 
                          triangle_ends = NULL, col_inf = NULL, col_sup = NULL,
                          colNA = NULL, color_fun = clim.palette(),
                          filled.continents = FALSE, coast_color = NULL, 
                          coast_width = 1,
                          dots = NULL, dot_symbol = 4, dot_size = 0.8,
                          intlat = 10, 
                          drawleg = TRUE, subsampleg = NULL, 
                          bar_extra_labels = NULL, draw_bar_ticks = TRUE, 
                          draw_separators = FALSE, triangle_ends_scale = 1,
                          bar_label_digits = 4, bar_label_scale = 1, 
                          units_scale = 1, bar_tick_scale = 1,
                          bar_extra_margin = rep(0, 4), 
                          boxlim = NULL, boxcol = "purple2", boxlwd = 5, 
                          margin_scale = rep(1, 4), title_scale = 1,
                          numbfig = NULL, fileout = NULL, 
                          width = 6, height = 5, size_units = 'in', 
                          res = 100, ...) {
  # Process the user graphical parameters that may be passed in the call
  ## Graphical parameters to exclude
  excludedArgs <- c("cex", "cex.main", "col", "fin", "lab", "las", "lwd", "mai", "mar", "mgp", "new", "pch", "ps")
  userArgs <- .FilterUserGraphicArgs(excludedArgs, ...)

  # If there is any filenames to store the graphics, process them
  # to select the right device 
  if (!is.null(fileout)) {
    deviceInfo <- .SelectDevice(fileout = fileout, width = width, height = height, units = size_units, res = res)
    saveToFile <- deviceInfo$fun
    fileout <- deviceInfo$files
  }

  # Preliminar check of dots, lon, lat
  if (!is.null(dots)) {
    if (!is.array(dots) || !(length(dim(dots)) %in% c(2, 3))) {
      stop("Parameter 'dots' must be a logical array with two or three dimensions.")
    }
    if (length(dim(dots)) == 2) {
      dim(dots) <- c(1, dim(dots))
    }
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

  # Check latlims
  if (!is.numeric(latlims) || length(latlims) != 2) {
    stop("Parameter 'latlims' must be a numeric vector with two elements.")
  }
  latlims <- sort(latlims)
  center_at <- 90 * sign(latlims[which.max(abs(latlims))])
  if (max(abs(latlims - center_at)) > 90 + 20) {
    stop("The range specified in 'latlims' is too wide. 110 degrees supported maximum.")
  }
  dlon <- median(lon[2:dims[1]] - lon[1:(dims[1] - 1)]) / 2
  dlat <- median(lat[2:dims[2]] - lat[1:(dims[2] - 1)]) / 2
  original_last_lat <- latlims[which.min(abs(latlims))]
  last_lat <- lat[which.min(abs(lat - original_last_lat))] - dlat * sign(center_at)
  latlims[which.min(abs(latlims))] <- last_lat

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

  # Check: brks, cols, subsampleg, bar_limits, color_fun, bar_extra_labels, draw_bar_ticks
  #        draw_separators, triangle_ends_scale, label_scale, units, units_scale, 
  #        bar_label_digits
  # Build: brks, cols, bar_limits, col_inf, col_sup
  var_limits <- c(min(var, na.rm = TRUE), max(var, na.rm = TRUE))
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

  # Check filled.continents
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

  # Check intlat
  if (!is.numeric(intlat)) {
    stop("Parameter 'intlat' must be numeric.")
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

  # Check numbfig
  if (!is.null(numbfig)) {
    if (!is.numeric(numbfig)) {
      stop("Parameter 'numbfig' must be numeric.")
    } else {
      numbfig <- round(numbfig)
      scale <- 1 / numbfig ** 0.3
      title_scale <- title_scale * scale
      margin_scale <- margin_scale * scale
      dot_size <- dot_size * scale
    }
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
  margins <- rep(0.2, 4) * margin_scale
  cex_title <- 2 * title_scale
  if (toptitle != '') {
    margins[3] <- margins[3] + cex_title + 1
  }
  bar_extra_margin[1] <- bar_extra_margin[1] + margins[1]
  bar_extra_margin[3] <- bar_extra_margin[3] + margins[3]

  if (drawleg) {
    layout(matrix(1:2, ncol = 2, nrow = 1), widths = c(8, 2))
  }
  # Load the user parameters
  par(userArgs)
  par(mar = margins, las = 0)
  coast <- map("world", interior = FALSE, projection = "stereographic", 
               orientation = c(center_at, 0, 0), fill = filled.continents,
               xlim = c(-180,180), ylim = latlims, wrap = TRUE, plot = FALSE)
  # Compute the bounding circle
  limit <- abs(mapproj::mapproject(0, last_lat, projection = 'stereographic',
                                   orientation = c(center_at, 0, 0))$y)
  for (i in 1:length(coast$x)) {
    distance <- sqrt(coast$x[i]**2 + coast$y[i]**2)
    if (!is.na(distance)) {
      if (distance > limit) {
        coast$x[i] <- coast$x[i] / distance * limit
        coast$y[i] <- coast$y[i] / distance * limit
      }
    }
  }
  xcircle <- c()
  ycircle <- c()
  for (i in 0:500) {
    xcircle <- c(xcircle, sin(2 * pi / 500 * i) * limit)
    ycircle <- c(ycircle, cos(2 * pi / 500 * i) * limit)
  }
  circle <- list(x = xcircle, y = ycircle)
  # Plot circle to set up device
  plot(circle, type= 'l', axes = FALSE, lwd = 1, col = gray(0.2), asp = 1,
       xlab = '', ylab = '', main = toptitle, cex.main = cex_title)
  col_inf_image <- ifelse(is.null(col_inf), colNA, col_inf)
  col_sup_image <- ifelse(is.null(col_sup), colNA, col_sup)
  # Draw the data polygons
  for (jx in 1:dims[1]) { 
    for (jy in 1:dims[2]) {
      if (lat[jy] >= latlims[1] && latlims[2] >= lat[jy]) {
        coord <- mapproj::mapproject(c(lon[jx] - dlon, lon[jx] + dlon,
                                       lon[jx] + dlon, lon[jx] - dlon),
                                     c(lat[jy] - dlat, lat[jy] - dlat,
                                       lat[jy] + dlat, lat[jy] + dlat))
        if (is.na(var[jx, jy] > 0)) {
          col <- colNA
        } else if (var[jx, jy] <= brks[1]) {
          col <- col_inf_image
        } else if (var[jx, jy] >= tail(brks, 1)) {
          col <- col_sup_image
        } else {
          ind <- which(brks[-1] >= var[jx, jy] & var[jx, jy] > brks[-length(brks)])
          col <- cols[ind]
        }
        polygon(coord, col = col, border = NA)
      }
    }
  }
  # Draw the dots
  if (!is.null(dots)) {
    numbfig <- 1  # for compatibility with PlotEquiMap code
    dots <- dots[, , which(lat >= latlims[1] & lat <= latlims[2]), drop = FALSE]
    data_avail <- !is.na(var[, which(lat >= latlims[1] & lat <= latlims[2]), drop = FALSE])
    for (counter in 1:(dim(dots)[1])) {
      points <- which(dots[counter, , ] & data_avail, arr.ind = TRUE)
      points_proj <- mapproj::mapproject(lon[points[, 1]], lat[points[, 2]])
      points(points_proj$x, points_proj$y,
             pch = dot_symbol[counter],
             cex = dot_size[counter] * 3 / sqrt(sqrt(sum(lat >= latlims[which.min(abs(latlims))]) * length(lon))),
             lwd = dot_size[counter] * 3 / sqrt(sqrt(sum(lat >= latlims[which.min(abs(latlims))]) * length(lon))))
    }
  }

  # Draw the continents, grid and bounding circle
  if (filled.continents) {
    old_lwd <- par('lwd')
    par(lwd = coast_width)
    polygon(coast, col = continent_color, border = coast_color)
    par(lwd = old_lwd)
  } else {
    lines(coast, col = coast_color, lwd = coast_width)
  }
  mapproj::map.grid(lim = c(-180, 180, latlims), nx = 18, 
                    ny = ceiling((latlims[2] - latlims[1]) / intlat),
                    col = 'lightgrey', labels = FALSE) 
  polygon(circle, border = 'black')
  # Draw boxes on the map
  if (!is.null(boxlim)) {
    counter <- 1
    for (box in boxlim) {
      if (box[1] > box[3]) {
        box[1] <- box[1] - 360
      }
      if (length(box) != 4) {
        stop(paste("The", counter, "st box defined in the parameter 'boxlim' is ill defined."))
      } else if (center_at == 90 && (box[2] < original_last_lat || 
                                     box[4] > center_at) ||
                 center_at == -90 && (box[4] > original_last_lat ||
                                      box[2] < center_at)) {
        stop(paste("The limits of the", counter, 
                   "st box defined in the parameter 'boxlim' are invalid."))
      } else {
        mapproj::map.grid(lim = c(box[1], box[3], box[2], box[4]), 
                          nx = 2, ny = 2, pretty = FALSE, 
                          col = boxcol[counter], lty = "solid",
                          lwd = boxlwd[counter], labels = FALSE)
      }
      counter <- counter + 1
    }
  }

  #
  #  Colorbar
  # ~~~~~~~~~~
  #
  if (drawleg) {
    ColorBar(brks, cols, TRUE, subsampleg, bar_limits, var_limits,
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
