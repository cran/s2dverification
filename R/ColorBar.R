ColorBar <- function(brks = NULL, cols = NULL, vertical = TRUE, 
                     subsampleg = NULL, bar_limits = NULL, var_limits = NULL, 
                     triangle_ends = NULL, col_inf = NULL, col_sup = NULL, 
                     color_fun = clim.palette(), plot = TRUE, 
                     draw_ticks = TRUE, draw_separators = FALSE, 
                     triangle_ends_scale = 1, extra_labels = NULL, 
                     title = NULL, title_scale = 1, 
                     label_scale = 1, tick_scale = 1, 
                     extra_margin = rep(0, 4), label_digits = 4, ...) {
  # Required checks
  if ((is.null(brks) || length(brks) < 2) && is.null(bar_limits) && is.null(var_limits)) {
    stop("At least one of 'brks' with the desired breaks, 'bar_limits' or ",
         "'var_limits' must be provided to generate the colour bar.")
  }

  # Check brks
  if (!is.null(brks)) {
    if (!is.numeric(brks)) {
      stop("Parameter 'brks' must be numeric if specified.")
    } else if (length(brks) > 1) {
      reorder <- sort(brks, index.return = TRUE)
      if (!is.null(cols)) {
        cols <- cols[reorder$ix[which(reorder$ix <= length(cols))]]
      }
      brks <- reorder$x
    } 
  }

  # Check bar_limits
  if (!is.null(bar_limits)) {
    if (!(all(is.na(bar_limits) | is.numeric(bar_limits)) && (length(bar_limits) == 2))) {
      stop("Parameter 'bar_limits' must be a vector of two numeric elements or NAs.")
    }
  }

  # Check var_limits
  if (!is.null(var_limits)) {
    if (!(is.numeric(var_limits) && (length(var_limits) == 2))) {
      stop("Parameter 'var_limits' must be a numeric vector of length 2.")
    } else if (any(is.na(var_limits))) {
      stop("Parameter 'var_limits' must not contain NA values.")
    } else if (any(is.infinite(var_limits))) {
      stop("Parameter 'var_limits' must not contain infinite values.")
    }
  }

  # Check cols
  if (!is.null(cols)) {
    if (!is.character(cols)) {
      stop("Parameter 'cols' must be a vector of character strings.")
    } else if (any(!sapply(cols, .IsColor))) {
      stop("Parameter 'cols' must contain valid colour identifiers.")
    }
  }

  # Check color_fun
  if (!is.function(color_fun)) {
    stop("Parameter 'color_fun' must be a colour-generator function.")
  }

  # Check integrity among brks, bar_limits and var_limits
  if (is.null(brks) || (length(brks) < 2)) {
    if (is.null(brks)) {
      if (is.null(cols)) {
        brks <- 21
      } else {
        brks <- length(cols) + 1
      }
    }
    if (is.null(bar_limits) || any(is.na(bar_limits))) {
      # var_limits is defined
      if (is.null(bar_limits)) {
        bar_limits <- c(NA, NA)
      }
      half_width <- 0.5 * (var_limits[2] - var_limits[1]) / (brks - 1)  
      bar_limits[which(is.na(bar_limits))] <- c(var_limits[1] - half_width, var_limits[2] + half_width)[which(is.na(bar_limits))]
      brks <- seq(bar_limits[1], bar_limits[2], length.out = brks)
    } else if (is.null(var_limits)) {
      # bar_limits is defined
      brks <- seq(bar_limits[1], bar_limits[2], length.out = brks)
      var_limits <- bar_limits
      var_limits[1] <- var_limits[1] + .Machine$double.xmin
    } else {
      # both bar_limits and var_limits are defined
      brks <- seq(bar_limits[1], bar_limits[2], length.out = brks)
    }
  } else if (is.null(bar_limits)) {
    if (is.null(var_limits)) {
      # brks is defined
      bar_limits <- c(head(brks, 1), tail(brks, 1))
      var_limits <- bar_limits
      var_limits[1] <- var_limits[1] + .Machine$double.xmin
    } else {
      # brks and var_limits are defined
      bar_limits <- c(head(brks, 1), tail(brks, 1))
    }
  } else {
    # brks and bar_limits are defined
    # or
    # brks, bar_limits and var_limits are defined
    if (head(brks, 1) != bar_limits[1] || tail(brks, 1) != bar_limits[2]) {
      stop("Parameters 'brks' and 'bar_limits' are inconsistent.")
    }
  }   
  
  # Check col_inf
  if (!is.null(col_inf)) {
    if (!.IsColor(col_inf)) {
      stop("Parameter 'col_inf' must be a valid colour identifier.")
    }
  }  

  # Check col_sup
  if (!is.null(col_sup)) {
    if (!.IsColor(col_sup)) {
      stop("Parameter 'col_sup' must be a valid colour identifier.")
    }
  }

  # Check triangle_ends
  if (!is.null(triangle_ends) && (!is.logical(triangle_ends) || length(triangle_ends) != 2)) {
    stop("Parameter 'triangle_ends' must be a logical vector with two elements.")
  }
  teflc <- triangle_ends_from_limit_cols <- c(!is.null(col_inf), !is.null(col_sup))
  if (is.null(triangle_ends) && is.null(col_inf) && is.null(col_sup)) {
    triangle_ends <- c(FALSE, FALSE)
    if (bar_limits[1] >= var_limits[1]) {
      triangle_ends[1] <- TRUE
    }
    if (bar_limits[2] < var_limits[2]) {
      triangle_ends[2] <- TRUE
    }
  } else if (!is.null(triangle_ends) && is.null(col_inf) && is.null(col_sup)) {
    triangle_ends <- triangle_ends
  } else if (is.null(triangle_ends) && (!is.null(col_inf) || !is.null(col_sup))) {
    triangle_ends <- teflc
  } else if (any(teflc != triangle_ends)) {
    if (!is.null(brks) && length(brks) > 1 && !is.null(cols) && length(cols) >= length(brks)) {
      triangle_ends <- teflc
    } else if (!is.null(cols)) {
      triangle_ends <- teflc
    } else {
      triangle_ends <- triangle_ends
    }
  }
  if (plot) {
    if ((bar_limits[1] >= var_limits[1]) && !triangle_ends[1]) {
      .warning("There are variable values smaller or equal to the lower limit ",
               "of the colour bar and the lower triangle end has been ",
               "disabled. These will be painted in the colour for NA values.")
    }
    if ((bar_limits[2] < var_limits[2]) && !triangle_ends[2]) {
      .warning("There are variable values greater than the higher limit ",
               "of the colour bar and the higher triangle end has been ",
               "disabled. These will be painted in the colour for NA values.")
    }
  }

  # Generate colours if needed
  if (is.null(cols)) {
    cols <- color_fun(length(brks) - 1 + sum(triangle_ends))
    attr_bk <- attributes(cols)
    if (triangle_ends[1]) {
      if (is.null(col_inf)) col_inf <- head(cols, 1)
      cols <- cols[-1]
    }
    if (triangle_ends[2]) {
      if (is.null(col_sup)) col_sup <- tail(cols, 1)
      cols <- cols[-length(cols)]
    }
    attributes(cols) <- attr_bk
  } else if ((length(cols) != (length(brks) - 1))) {
    stop("Incorrect number of 'brks' and 'cols'. There must be one more break than the number of colours.")
  }

  # Check vertical
  if (!is.logical(vertical)) {
    stop("Parameter 'vertical' must be TRUE or FALSE.")
  }

  # Check extra_labels
  if (is.null(extra_labels)) {
    extra_labels <- numeric(0)
  }
  if (!is.numeric(extra_labels)) {
    stop("Parameter 'extra_labels' must be numeric.")
  } else {
    if (any(extra_labels > bar_limits[2]) || any(extra_labels < bar_limits[1])) {
      stop("Parameter 'extra_labels' must not contain ticks beyond the color bar limits.")
    }
  }
  extra_labels <- sort(extra_labels)

  # Check subsampleg
  primes <- function(x) {
    # Courtesy of Chase. See http://stackoverflow.com/questions/6424856/r-function-for-returning-all-factors
    x <- as.integer(x)
    div <- seq_len(abs(x))
    factors <- div[x %% div == 0L]
    factors <- list(neg = -factors, pos = factors)
    return(factors)
  }
  remove_final_tick <- FALSE
  added_final_tick <- TRUE
  if (is.null(subsampleg)) {
    subsampleg <- 1
    while (length(brks) / subsampleg > 15 - 1) {
      next_factor <- primes((length(brks) - 1) / subsampleg)$pos
      next_factor <- next_factor[length(next_factor) - ifelse(length(next_factor) > 2, 1, 0)]
      subsampleg <- subsampleg * next_factor
    }
    if (subsampleg > (length(brks) - 1) / 4) {
      subsampleg <- max(1, round(length(brks) / 4))
      extra_labels <- c(extra_labels, bar_limits[2])
      added_final_tick <- TRUE
      if ((length(brks) - 1) %% subsampleg < (length(brks) - 1) / 4 / 2) {
        remove_final_tick <- TRUE
      }
    }
  } else if (!is.numeric(subsampleg)) {
    stop("Parameter 'subsampleg' must be numeric.")
  }
  subsampleg <- round(subsampleg)

  # Check plot
  if (!is.logical(plot)) {
    stop("Parameter 'plot' must be logical.")
  }

  # Check draw_separators
  if (!is.logical(draw_separators)) {
    stop("Parameter 'draw_separators' must be logical.")
  }

  # Check triangle_ends_scale
  if (!is.numeric(triangle_ends_scale)) {
    stop("Parameter 'triangle_ends_scale' must be numeric.")
  }

  # Check draw_ticks
  if (!is.logical(draw_ticks)) {
    stop("Parameter 'draw_ticks' must be logical.")
  }

  # Check title
  if (is.null(title)) {
    title <- ''
  }
  if (!is.character(title)) {
    stop("Parameter 'title' must be a character string.")
  }

  # Check title_scale
  if (!is.numeric(title_scale)) {
    stop("Parameter 'title_scale' must be numeric.")
  }

  # Check label_scale
  if (!is.numeric(label_scale)) {
    stop("Parameter 'label_scale' must be numeric.")
  }

  # Check tick_scale
  if (!is.numeric(tick_scale)) {
    stop("Parameter 'tick_scale' must be numeric.")
  }

  # Check extra_margin
  if (!is.numeric(extra_margin) || length(extra_margin) != 4) {
    stop("Parameter 'extra_margin' must be a numeric vector of length 4.")
  }

  # Check label_digits
  if (!is.numeric(label_digits)) {
    stop("Parameter 'label_digits' must be numeric.")
  }
  label_digits <- round(label_digits)

  # Process the user graphical parameters that may be passed in the call
  ## Graphical parameters to exclude
  excludedArgs <- c("cex", "cex.axis", "col", "lab", "las", "mar", "mgp", "new", "ps")
  userArgs <- .FilterUserGraphicArgs(excludedArgs, ...)

  #
  #  Plotting colorbar
  # ~~~~~~~~~~~~~~~~~~~
  #
  if (plot) {
    pars_to_save <- c('mar', 'cex', names(userArgs), 'mai', 'mgp', 'las', 'xpd')
    saved_pars <- par(pars_to_save) 
    par(mar = c(0, 0, 0, 0), cex = 1)
    image(1, 1, t(t(1)), col = rgb(0, 0, 0, 0), axes = FALSE, xlab = '', ylab = '')
    # Get the availale space
    figure_size <- par('fin')
    cs <- par('csi')
    # This allows us to assume we always want to plot horizontally
    if (vertical) {
      figure_size <- rev(figure_size)
    }
#    pannel_to_redraw <- par('mfg')
#    .SwitchToFigure(pannel_to_redraw[1], pannel_to_redraw[2])
    # Load the user parameters
    par(new = TRUE)
    par(userArgs)
    # Set up color bar plot region
    margins <- c(0.0, 0, 0.0, 0)
    cex_title <- 1 * title_scale
    cex_labels <- 0.9 * label_scale
    cex_ticks <- -0.3 * tick_scale
    spaceticklab <- max(-cex_ticks, 0)
    if (vertical) {
      margins[1] <- margins[1] + (1.2 * cex_labels * 3 + spaceticklab) * cs
      margins <- margins + extra_margin[c(4, 1:3)] * cs
    } else {
      margins[1] <- margins[1] + (1.2 * cex_labels * 1 + spaceticklab) * cs
      margins <- margins + extra_margin * cs
    }
    if (title != '') {
      margins[3] <- margins[3] + (1.0 * cex_title) * cs
    }
    margins[3] <- margins[3] + sqrt(figure_size[2] / (margins[1] + margins[3])) * 
                               figure_size[2] / 6 * ifelse(title != '', 0.5, 0.8)
    # Set side margins
    margins[2] <- margins[2] + figure_size[1] / 16
    margins[4] <- margins[4] + figure_size[1] / 16
    triangle_ends_prop <- 1 / 32 * triangle_ends_scale 
    triangle_ends_cex <- triangle_ends_prop * figure_size[2]
    if (triangle_ends[1]) {
      margins[2] <- margins[2] + triangle_ends_cex
    }
    if (triangle_ends[2]) {
      margins[4] <- margins[4] + triangle_ends_cex
    }
    ncols <- length(cols)
    # Set up the points of triangles
    # Compute the proportion of horiz. space occupied by one plot unit
    prop_unit <- (1 - (margins[2] + margins[4]) / figure_size[1]) / ncols
    # Convert triangle height to plot inits
    triangle_height <- triangle_ends_prop / prop_unit
    left_triangle <- list(x = c(1, 1 - triangle_height, 1) - 0.5,
                          y = c(1.4, 1, 0.6))
    right_triangle <- list(x = c(ncols, ncols + triangle_height, ncols) + 0.5,
                           y = c(1.4, 1, 0.6))
    # Draw the color squares and title
    if (vertical) {
      par(mai = c(margins[2:4], margins[1]), 
          mgp = c(0, spaceticklab + 0.2, 0), las = 1)
      d <- 4
      image(1, 1:ncols, t(1:ncols), axes = FALSE, col = cols, 
            xlab = '', ylab = '')
      title(ylab = title, line = cex_title * (0.2 + 0.1), cex.lab = cex_title)
      # Draw top and bottom border lines    
      lines(c(0.6, 0.6), c(1 - 0.5, ncols + 0.5))
      lines(c(1.4, 1.4), c(1 - 0.5, ncols + 0.5))
      # Rotate triangles
      names(left_triangle) <- rev(names(left_triangle))
      names(right_triangle) <- rev(names(right_triangle))
    } else {
      # The term - cex_labels / 4 * (3 / cex_labels - 1) was found by
      # try and error
      par(mai = margins, 
          mgp = c(0, cex_labels / 2 + spaceticklab 
                     - cex_labels / 4 * (3 / cex_labels - 1), 0),
          las = 1)
      d <- 1
      image(1:ncols, 1, t(t(1:ncols)), axes = FALSE, col = cols, 
            xlab = '', ylab = '')
      title(title, line = cex_title * (0.2 + 0.1), cex.main = cex_title)
      # Draw top and bottom border lines    
      lines(c(1 - 0.5, ncols + 0.5), c(0.6, 0.6))
      lines(c(1 - 0.5, ncols + 0.5), c(1.4, 1.4))
      tick_length <- -0.4
    }
    # Draw the triangles
    par(xpd = TRUE)
    if (triangle_ends[1]) {
      # Draw left triangle
      polygon(left_triangle$x, left_triangle$y, col = col_inf, border = NA)
      lines(left_triangle$x, left_triangle$y)       
    }
    if (triangle_ends[2]) {
      # Draw right triangle
      polygon(right_triangle$x, right_triangle$y, col = col_sup, border = NA)
      lines(right_triangle$x, right_triangle$y)
    }
    par(xpd = FALSE)

    # Put the separators
    if (vertical) {
      if (draw_separators) {
        for (i in 1:(ncols - 1)) {
          lines(c(0.6, 1.4), c(i, i) + 0.5)
        }
      }
      if (draw_separators || is.null(col_inf)) {
        lines(c(0.6, 1.4), c(0.5, 0.5))
      }
      if (draw_separators || is.null(col_sup)) {
        lines(c(0.6, 1.4), c(ncols + 0.5, ncols + 0.5))
      }
    } else {
      if (draw_separators) {
        for (i in 1:(ncols - 1)) {
          lines(c(i, i) + 0.5, c(0.6, 1.4))
        }
      }
      if (draw_separators || is.null(col_inf)) {
        lines(c(0.5, 0.5), c(0.6, 1.4))
      }
      if (draw_separators || is.null(col_sup)) {
        lines(c(ncols + 0.5, ncols + 0.5), c(0.6, 1.4))
      }
    }
    # Put the ticks
    plot_range <- length(brks) - 1
    var_range <- tail(brks, 1) - head(brks, 1)
    extra_labels_at <- ((extra_labels - head(brks, 1)) / var_range) * plot_range + 0.5
    at <- seq(1, length(brks), subsampleg)
    labels <- brks[at]
    # Getting rid of next-to-last tick if too close to last one
    if (remove_final_tick) {
      at <- at[-length(at)]
      labels <- labels[-length(labels)]
    }
    labels <- signif(labels, label_digits)
    if (added_final_tick) {
      extra_labels[length(extra_labels)] <- signif(tail(extra_labels, 1), label_digits)
    }
    at <- at - 0.5
    at <- c(at, extra_labels_at)
    labels <- c(labels, extra_labels)
    tick_reorder <- sort(at, index.return = TRUE)
    at <- tick_reorder$x
    labels <- labels[tick_reorder$ix]
    axis(d, at = at, tick = draw_ticks, labels = labels, cex.axis = cex_labels, tcl = cex_ticks)
    par(saved_pars)
  }
  invisible(list(brks = brks, cols = cols, col_inf = col_inf, col_sup = col_sup))
}
