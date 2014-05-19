ColorBar <- function(brks, cols = NULL, vert = TRUE, subsampleg = 1) {
  # Creates a horizontal or vertical colorbar to introduce in multipanels.
  #
  # Args:
  #   brks: Levels.
  #   cols: List of colours, optional.
  #   vert: TRUE/FALSE for vertical/horizontal colorbar.
  #   kharin: Supsampling factor of the interval between ticks on colorbar.
  #           Default: 1 = every level
  #
  # Returns:
  #   This function returns nothing
  #
  # History:
  #   1.0  #  2012-04  (V. Guemas, vguemas@ic3.cat)  #  Original code

  #
  #
  #  Input arguments
  # ~~~~~~~~~~~~~~~~~
  #
  if (is.null(cols) == TRUE) {
    nlev <- length(brks) - 1
    cols <- rainbow(nlev)
  } else {
    if (length(cols) != (length(brks) - 1)) {
      stop("Inconsistent colour levels / list of colours")
    }
  }
  
  # 
  #  Plotting colorbar
  # ~~~~~~~~~~~~~~~~~~~
  #
  if (vert) {
    par(mar = c(1, 1, 1, 2.5), mgp = c(1, 1, 0), las = 1, cex = 1.2)
    image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
          xlab = '', ylab = '')
    box()
    axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE, 
         labels = brks[seq(1, length(brks), subsampleg)])
  } else {
    par(mar = c(1.5, 1, 1, 1), mgp = c(1.5, 0.3, 0), las = 1, cex = 1.2)
    image(1:length(cols), 1, t(t(1:length(cols))), axes = FALSE, col = cols,
          xlab = '', ylab = '')
    box()
    axis(1, at = seq(0.5, length(brks) - 0.5, subsampleg), 
         labels = brks[seq(1, length(brks), subsampleg)])
  }
}
