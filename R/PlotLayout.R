PlotLayout <- function(fun, plot_dims, var, ..., special_args = NULL, 
                       nrow = NULL, ncol = NULL, toptitle = NULL, 
                       row_titles = NULL, col_titles = NULL, bar_scale = 1, 
                       title_scale = 1, title_margin_scale = 1,
                       title_left_shift_scale = 1,
                       subtitle_scale = 1, subtitle_margin_scale = 1,
                       brks = NULL, cols = NULL, drawleg = 'S', titles = NULL, 
                       subsampleg = NULL, bar_limits = NULL, 
                       triangle_ends = NULL, col_inf = NULL, col_sup = NULL,
                       color_fun = clim.colors, 
                       draw_bar_ticks = TRUE, draw_separators = FALSE, 
                       triangle_ends_scale = 1, bar_extra_labels = NULL, 
                       units = NULL, units_scale = 1, bar_label_scale = 1, 
                       bar_tick_scale = 1, bar_extra_margin = rep(0, 4), 
                       bar_left_shift_scale = 1, bar_label_digits = 4, 
                       extra_margin = rep(0, 4),
                       fileout = NULL, width = NULL, height = NULL, 
                       size_units = 'in', res = 100, close_device = TRUE) {
  # If there is any filenames to store the graphics, process them
  # to select the right device 
  if (!is.null(fileout)) {
    deviceInfo <- .SelectDevice(fileout = fileout, width = width, height = height, units = size_units, res = res)
    saveToFile <- deviceInfo$fun
    fileout <- deviceInfo$files
  }

  is_single_na <- function(x) ifelse(length(x) > 1, FALSE, is.na(x))
  # Check var
  if (is.array(var) || (is_single_na(var))) {
    var <- list(var)
  } else if (is.list(var)) {
    if (!all(sapply(var, is.array) | sapply(var, is_single_na))) {
      stop("Parameter 'var' must be an array or a list of arrays (or NA values).")
    }
  } else {
    stop("Parameter 'var' must be an array or a list of arrays.")
  }

  # Check fun
  if (length(fun) == 1) {
    if (is.function(fun)) {
      fun <- as.character(substitute(fun))
    }
    if (is.character(fun)) {
      fun <- rep(fun, length(var))
    }
  }
  if (!is.character(fun) || (length(fun) != length(var))) {
    stop("Parameter 'fun' must be a single function or a vector of function names, one for each array provided in parameter 'var'.")
  }

  # Check special_args
  if (!is.null(special_args)) {
    if (!is.list(special_args) || any(!sapply(special_args, is.list))) {
      stop("Parameter 'special_args' must be a list of lists.")
    } else if (length(special_args) != length(var)) {
      stop("Parameter 'special_args' must contain a list of special arguments for each array provided in 'var'.")
    }
  }

  # Check plot_dims
  if (is.character(plot_dims) || is.numeric(plot_dims)) {
    plot_dims <- replicate(length(var), plot_dims, simplify = FALSE)
  }
  if (!is.list(plot_dims) || !all(sapply(plot_dims, is.character) | sapply(plot_dims, is.numeric)) ||
      (length(plot_dims) != length(var))) {
    stop("Parameter 'plot_dims' must contain a single numeric or character vector with dimension identifiers or a vector for each array provided in parameter 'var'.")
  }

  # Check nrow
  if (!is.null(nrow)) {
    if (!is.numeric(nrow)) {
      stop("Parameter 'nrow' must be numeric or NULL.")
    }
    nrow <- round(nrow)
  }

  # Check ncol
  if (!is.null(ncol)) {
    if (!is.numeric(ncol)) {
      stop("Parameter 'ncol' must be numeric or NULL.")
    }
    ncol <- round(ncol)
  }

  # Check toptitle
  if (is.null(toptitle) || is.na(toptitle)) {
    toptitle <- ''
  }
  if (!is.character(toptitle)) {
    stop("Parameter 'toptitle' must be a character string.")
  }

  # Check row_titles
  if (!is.null(row_titles)) {
    if (!is.character(row_titles)) {
      stop("Parameter 'row_titles' must be a vector of character strings.")
    }
  }

  # Check col_titles
  if (!is.null(row_titles)) {
    if (!is.character(row_titles)) {
      stop("Parameter 'row_titles' must be a vector of character strings.")
    }
  }

  # Check drawleg
  if (is.character(drawleg)) {
    if (drawleg %in% c('up', 'u', 'U', 'top', 't', 'T', 'north', 'n', 'N')) {
      drawleg <- 'N'
    } else if (drawleg %in% c('down', 'd', 'D', 'bottom', 'b', 'B', 'south', 's', 'S')) {
      drawleg <- 'S'
    } else if (drawleg %in% c('right', 'r', 'R', 'east', 'e', 'E')) {
      drawleg <- 'E'
    } else if (drawleg %in% c('left', 'l', 'L', 'west', 'w', 'W')) {
      drawleg <- 'W'
    } else {
      stop("Parameter 'drawleg' must be either TRUE, FALSE or a valid identifier of a position (see ?PlotMultiMap).")
    }
  } else if (!is.logical(drawleg)) {
    stop("Parameter 'drawleg' must be either TRUE, FALSE or a valid identifier of a position (see ?PlotMultiMap).")
  }
  if (drawleg != FALSE && all(sapply(var, is_single_na)) && 
      (is.null(brks) || length(brks) < 2)) {
    stop("Either data arrays in 'var' or breaks in 'brks' must be provided if 'drawleg' is requested.")
  }

  # Check the rest of parameters (unless the user simply wants to build an empty layout)
  var_limits <- NULL
  if (!all(sapply(var, is_single_na))) {
    var_limits <- c(min(unlist(var), na.rm = TRUE), max(unlist(var), na.rm = TRUE))
    if ((any(is.infinite(var_limits)) || var_limits[1] == var_limits[2]))  {
      stop("Arrays in parameter 'var' must contain at least 2 different values.")
    }
  }
  colorbar <- ColorBar(brks, cols, FALSE, subsampleg, bar_limits,
                       var_limits, triangle_ends, col_inf, col_sup, color_fun, 
                       plot = FALSE, draw_bar_ticks, 
                       draw_separators, triangle_ends_scale, bar_extra_labels,
                       units, units_scale, bar_label_scale, bar_tick_scale,
                       bar_extra_margin, bar_label_digits)
 
  # Check bar_scale
  if (!is.numeric(bar_scale)) {
    stop("Parameter 'bar_scale' must be numeric.")
  }

  # Check bar_left_shift_scale
  if (!is.numeric(bar_left_shift_scale)) {
    stop("Parameter 'bar_left_shift_scale' must be numeric.")
  }

  # Check title_scale
  if (!is.numeric(title_scale)) {
    stop("Parameter 'title_scale' must be numeric.")
  }

  # Check title_margin_scale
  if (!is.numeric(title_margin_scale)) {
    stop("Parameter 'title_margin_scale' must be numeric.")
  }

  # Check title_left_shift_scale
  if (!is.numeric(title_left_shift_scale)) {
    stop("Parameter 'title_left_shift_scale' must be numeric.")
  }

  # Check subtitle_scale
  if (!is.numeric(subtitle_scale)) {
    stop("Parameter 'subtite_scale' must be numeric.")
  }

  # Check subtitle_margin_scale
  if (!is.numeric(subtitle_margin_scale)) {
    stop("Parameter 'subtite_margin_scale' must be numeric.")
  }

  # Check titles
  if (!all(sapply(titles, is.character))) {
    stop("Parameter 'titles' must be a vector of character strings.")
  }

  # Check extra_margin
  if (!is.numeric(extra_margin) || length(extra_margin) != 4) {
    stop("Parameter 'extra_margin' must be a numeric vector with 4 elements.")
  }

  # Check width
  if (is.null(width)) {
    if (is.null(fileout)) {
      width <- 7
    } else {
      width <- 11
    }
  }
  if (!is.numeric(width)) {
    stop("Parameter 'width' must be numeric.")
  }

  # Check height
  if (is.null(height)) {
    if (is.null(fileout)) {
      height <- 7
    } else {
      height <- 8
    }
  }
  if (!is.numeric(height)) {
    stop("Parameter 'height' must be numeric.")
  }

  # Check close_device
  if (!is.logical(close_device)) {
    stop("Parameter 'close_device' must be logical.")
  }

  # Count the total number of maps and reorder each array of maps to have the lat and lon dimensions at the end.
  n_plots <- 0
  plot_array_i <- 1
  for (plot_array in var) {
    if (is_single_na(plot_array)) {
      n_plots <- n_plots + 1
    } else {
      dim_ids <- plot_dims[[plot_array_i]]
      if (is.character(dim_ids)) {
        dimnames <- NULL
        if (!is.null(names(dim(plot_array)))) {
          dimnames <- names(dim(plot_array))
        } else if (!is.null(attr(plot_array, 'dimensions'))) {
          dimnames <- attr(plot_array, 'dimensions')
        }
        if (!is.null(dimnames)) {
          if (any(!sapply(dim_ids, `%in%`, dimnames))) {
            stop("All arrays provided in parameter 'var' must have all the dimensions in 'plot_dims'.")
          }
          dim_ids <- sapply(dim_ids, function(x) which(dimnames == x)[1])
          var[[plot_array_i]] <- .aperm2(var[[plot_array_i]], c((1:length(dim(plot_array)))[-dim_ids], dim_ids))
        } else {
          .warning(paste0("Assuming the ", plot_array_i, "th array provided in 'var' has 'plot_dims' as last dimensions (right-most)."))
          dims <- tail(c(rep(1, length(dim_ids)), dim(plot_array)), length(dim_ids))
          dim_ids <- tail(1:length(dim(plot_array)), length(dim_ids))
          if (length(dim(var[[plot_array_i]])) < length(dims)) {
            dim(var[[plot_array_i]]) <- dims
          }
        }
      } else if (any(dim_ids > length(dim(plot_array)))) {
        stop("Parameter 'plot_dims' contains dimension identifiers out of range.")
      }
      n_plots <- n_plots + prod(dim(plot_array)[-dim_ids])
      #n_plots <- n_plots + prod(head(c(rep(1, length(dim_ids)), dim(plot_array)), length(dim(plot_array))))
      if (length(dim(var[[plot_array_i]])) == length(dim_ids)) {
        dim(var[[plot_array_i]]) <- c(1, dim(var[[plot_array_i]]))
        dim_ids <- dim_ids + 1
      }
      plot_dims[[plot_array_i]] <- dim_ids
    }
    plot_array_i <- plot_array_i + 1
  }
  if (is.null(nrow) && is.null(ncol)) {
    ncol <- ceiling(sqrt(n_plots))
    nrow <- ceiling(n_plots/ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(n_plots/nrow)
  } else if (is.null(nrow)) {
    nrow <- ceiling(n_plots/ncol)
  } else if (nrow * ncol < n_plots) {
    stop("There are more arrays to plot in 'var' than cells defined by 'nrow' x 'ncol'.")
  }

  if (is.logical(drawleg) && drawleg) {
    if (nrow > ncol) {
      drawleg <- 'S'
    } else {
      drawleg <- 'E'
    }
  }
  vertical <- drawleg %in% c('E', 'W')

  # Open connection to graphical device
  if (!is.null(fileout)) {
    saveToFile(fileout)
  } else if (names(dev.cur()) == 'null device') {
    dev.new(units = size_units, res = res, width = width, height = height)
  } else if (prod(par('mfrow')) > 1) {
    dev.new(units = units, res = res, width = width, height = height)
  }

  # Take size of device and set up layout:
  # ---------------------------------------------
  # |0000000000000000000000000000000000000000000|
  # |0000000000000000 TOP TITLE 0000000000000000|
  # |0000000000000000000000000000000000000000000|
  # |-------------------------------------------|
  # |00000|0000000000000000000000000000000000000|
  # |00000|000000000000 ROW TITLES 0000000000000|
  # |00000|0000000000000000000000000000000000000|
  # |00000|-------------------------------------|
  # |0   0|222222222222222222|333333333333333333|
  # |0 C 0|222222222222222222|333333333333333333|
  # |0 O 0|222222222222222222|333333333333333333|
  # |0 L 0|2222 FIGURE 1 2222|3333 FIGURE 2 3333|
  # |0   0|222222222222222222|333333333333333333|
  # |0 T 0|222222222222222222|333333333333333333|
  # |0 I 0|222222222222222222|333333333333333333|
  # |0 T 0|-------------------------------------|
  # |0 L 0|444444444444444444|555555555555555555|
  # |0 S 0|444444444444444444|555555555555555555|
  # |0   0|444444444444444444|555555555555555555|
  # |00000|4444 FIGURE 3 4444|5555 FIGURE 4 5555|
  # |00000|444444444444444444|555555555555555555|
  # |00000|444444444444444444|555555555555555555|
  # |00000|444444444444444444|555555555555555555|
  # |-------------------------------------------|
  # |1111111111111111111111111111111111111111111|
  # |1111111111111111 COLOR BAR 1111111111111111|
  # |1111111111111111111111111111111111111111111|
  # ---------------------------------------------
  device_size <- par('din')
  device_size[1] <- device_size[1] - sum(extra_margin[c(2, 4)])
  device_size[2] <- device_size[2] - sum(extra_margin[c(1, 3)])
  cs <- char_size <- par('csi')
  title_cex <- 2.5 * title_scale
  title_margin <- 0.5 * title_cex * title_margin_scale
  subtitle_cex <- 1.5 * subtitle_scale
  subtitle_margin <- 0.5 * sqrt(nrow * ncol) * subtitle_cex * subtitle_margin_scale
  mat_layout <- 1:(nrow * ncol) + ifelse(drawleg != FALSE, 1, 0)
  mat_layout <- matrix(mat_layout, nrow, ncol, byrow = TRUE)
  fsu <- figure_size_units <- 10  # unitless
  widths <- rep(fsu, ncol)
  heights <- rep(fsu, nrow)
  n_figures <- nrow * ncol
  if (length(row_titles) > 0) {
    mat_layout <- cbind(rep(0, dim(mat_layout)[1]), mat_layout)
    widths <- c(((subtitle_cex + subtitle_margin / 2) * cs / device_size[1]) * ncol * fsu, widths)
  }
  if (length(col_titles) > 0) {
    mat_layout <- rbind(rep(0, dim(mat_layout)[2]), mat_layout)
    heights <- c(((subtitle_cex + subtitle_margin) * cs / device_size[2]) * nrow * fsu, heights)
  }
  if (drawleg != FALSE) {
    if (drawleg == 'N') {
      mat_layout <- rbind(rep(1, dim(mat_layout)[2]), mat_layout)
      heights <- c(round(bar_scale * 2 * nrow), heights)
    } else if (drawleg == 'S') {
      mat_layout <- rbind(mat_layout, rep(1, dim(mat_layout)[2]))
      heights <- c(heights, round(bar_scale * 2 * nrow))
    } else if (drawleg == 'W') {
      mat_layout <- cbind(rep(1, dim(mat_layout)[1]), mat_layout)
      widths <- c(round(bar_scale * 3 * ncol), widths)
    } else if (drawleg == 'E') {
      mat_layout <- cbind(mat_layout, rep(1, dim(mat_layout)[1]))
      widths <- c(widths, round(bar_scale * 3 * ncol))
    }
    n_figures <- n_figures + 1
  }
  if (toptitle != '') {
    mat_layout <- rbind(rep(0, dim(mat_layout)[2]), mat_layout)
    heights <-  c(((title_cex + title_margin) * cs / device_size[2]) * nrow * fsu, heights)
  }
  par(oma = extra_margin)
  layout(mat_layout, widths, heights)
  # Draw the color bar
  if (drawleg != FALSE) {
    if (length(row_titles) > 0) {
      bar_extra_margin[2] <- bar_extra_margin[2] + (subtitle_cex + subtitle_margin / 2) * 
                                                   bar_left_shift_scale
    }
    ColorBar(colorbar$brks, colorbar$cols, vertical, subsampleg, 
             bar_limits, var_limits,
             triangle_ends = triangle_ends, col_inf = colorbar$col_inf,
             col_sup = colorbar$col_sup, color_fun, plot = TRUE, draw_bar_ticks, 
             draw_separators, triangle_ends_scale, bar_extra_labels,
             units, units_scale, bar_label_scale, bar_tick_scale, 
             bar_extra_margin, bar_label_digits)
  }

  # Draw titles
  if (toptitle != '' || length(col_titles) > 0 || length(row_titles) > 0) {
    plot(0, type = 'n', ann = FALSE, axes = FALSE, xaxs = 'i', yaxs = 'i', 
         xlim = c(0, 1), ylim = c(0, 1))
    width_lines <- par('fin')[1] / par('csi')
    plot_lines <- par('pin')[1] / par('csi')
    plot_range <- par('xaxp')[2] - par('xaxp')[1]
    size_units_per_line <- plot_range / plot_lines
    if (toptitle != '') {
      title_x_center <- par('xaxp')[1] - par('mar')[2] * size_units_per_line + 
                        ncol * width_lines * size_units_per_line / 2
      if (length(row_titles) > 0) {
        title_x_center <- title_x_center - (1 - title_left_shift_scale) * 
                          (subtitle_cex + subtitle_margin) / 2 * size_units_per_line
      }
      title_y_center <- par('mar')[3] + (title_margin + title_cex) / 2
      if (length(col_titles > 0)) {
        title_y_center <- title_y_center + (subtitle_margin + subtitle_cex)
      }
      mtext(toptitle, cex = title_cex, line = title_y_center, at = title_x_center,
            padj = 0.5)
    }
    if (length(col_titles) > 0) {
      t_x_center <- par('xaxp')[1] - par('mar')[2] * size_units_per_line
      for (t in 1:ncol) {
        mtext(col_titles[t], cex = subtitle_cex,
              line = par('mar')[3] + (subtitle_margin + subtitle_cex) / 2,
              at = t_x_center + (t - 0.5) * width_lines * size_units_per_line,
              padj = 0.5)
      }
    }
    height_lines <- par('fin')[2] / par('csi')
    plot_lines <- par('pin')[2] / par('csi')
    plot_range <- par('yaxp')[2] - par('yaxp')[1]
    size_units_per_line <- plot_range / plot_lines
    if (length(row_titles) > 0) {
      t_y_center <- par('yaxp')[1] - par('mar')[1] * size_units_per_line
      for (t in 1:nrow) {
        mtext(row_titles[t], cex = subtitle_cex,
              line = par('mar')[2] + (subtitle_margin + subtitle_cex) / 2, 
              at = t_y_center - (t - 1.5) * height_lines * size_units_per_line, 
              padj = 0.5, side = 2)
      }
    }
    par(new = TRUE)
  }

  array_number <- 1
  plot_number <- 1
  # For each array provided in var
  lapply(var, function(x) {
    if (is_single_na(x)) {
      if (!all(sapply(var[array_number:length(var)], is_single_na))) {
        plot.new()
        par(new = FALSE)
      }
      plot_number <<- plot_number + 1
    } else {
      if (is.character(plot_dims[[array_number]])) {
        plot_dim_indices <- which(names(dim(x)) %in% plot_dims[[array_number]])
      } else {
        plot_dim_indices <- plot_dims[[array_number]]
      }
      # For each of the arrays provided in that array
      apply(x, (1:length(dim(x)))[-plot_dim_indices], 
            function(y) {
        # Do the plot
        fun_args <- c(list(y, toptitle = titles[plot_number]), list(...), 
                      special_args[[array_number]])
        funct <- fun[[array_number]]
        if (fun[[array_number]] %in% c('PlotEquiMap', 'PlotStereoMap', 'PlotSection')) {
          fun_args <- c(fun_args, list(brks = colorbar$brks, cols = colorbar$cols, 
                                       col_inf = colorbar$col_inf, 
                                       col_sup = colorbar$col_sup, 
                                       drawleg = FALSE)) 
        }
        do.call(fun[[array_number]], fun_args)
        plot_number <<- plot_number + 1
            })
    }
    array_number <<- array_number + 1
  })

  # If the graphic was saved to file, close the connection with the device
  if (!is.null(fileout) && close_device) dev.off()

  invisible(list(brks = colorbar$brks, cols = colorbar$cols, 
                 col_inf = colorbar$col_inf, col_sup = colorbar$col_sup,
                 layout_matrix = mat_layout))
}
