clim.colors <- function(n, palette = "bluered") {
  clim.palette(palette)(n)
}

clim.palette <- function(palette = "bluered") {
  if (palette == "bluered") {
    colorbar <- colorRampPalette(rev(c("#67001f", "#b2182b", "#d6604d",
                                       "#f4a582", "#fddbc7", "#f7f7f7", 
                                       "#d1e5f0", "#92c5de", "#4393c3", 
                                       "#2166ac", "#053061")))
    attr(colorbar, 'na_color') <- 'pink'
  } else if (palette == "redblue") {
    colorbar <- colorRampPalette(c("#67001f", "#b2182b", "#d6604d",
                                       "#f4a582", "#fddbc7", "#f7f7f7", 
                                       "#d1e5f0", "#92c5de", "#4393c3", 
                                       "#2166ac", "#053061"))
    attr(colorbar, 'na_color') <- 'pink'
  } else if (palette == "yellowred") {
    colorbar <- colorRampPalette(c("#ffffcc", "#ffeda0", "#fed976",
                                   "#feb24c", "#fd8d3c", "#fc4e2a",
                                   "#e31a1c", "#bd0026", "#800026"))
    attr(colorbar, 'na_color') <- 'pink'
  } else if (palette == "redyellow") {
    colorbar <- colorRampPalette(rev(c("#ffffcc", "#ffeda0", "#fed976",
                                   "#feb24c", "#fd8d3c", "#fc4e2a",
                                   "#e31a1c", "#bd0026", "#800026")))
    attr(colorbar, 'na_color') <- 'pink'
  } else {
    stop("Parameter 'palette' must be one of 'bluered', 'redblue', 'yellowred' or 'redyellow'.")
  }
  colorbar
}
