\name{ColorBar}
\alias{ColorBar}
\title{Draws Color Bar}
\description{
  Creates a horizontal or vertical colorbar to introduce in multipanels.
}
\usage{
  ColorBar(brks, cols = NULL, vert = TRUE, subsampleg = 1)
}
\arguments{
  \item{brks}{
    Levels.
  }
  \item{cols}{
    List of colours, optional.
  }
  \item{vert}{
    TRUE/FALSE for vertical/horizontal colorbar.
  }
  \item{subsampleg}{
    Supsampling factor of the interval between ticks on the colorbar.
    Default: 1 = every level
  }
}
\examples{
  cols <- c("dodgerblue4", "dodgerblue1", "forestgreen", "yellowgreen", "white",
            "white", "yellow", "orange", "red", "saddlebrown")
  lims <- seq(-1, 1, 0.2)
  ColorBar(lims, cols)
}
\author{
  History:
  0.1  -  2012-04  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  0.2  -  2013-04  (I. Andreu-Burillo, \email{isabel.andreu-burillo@ic3.cat})  -  vert option
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{dplot}
