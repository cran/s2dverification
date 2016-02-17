\name{clim.colors}
\alias{clim.colors}
\title{Generate Climate Color Bar}
\description{
Generates a color bar with color ranges useful in climate temperature variable plotting.\cr
The original colors are:\cr
  c("dodgerblue4", "dodgerblue1", "forestgreen", "yellowgreen", "white",\cr
    "white", "yellow", "orange", "red", "saddlebrown")
}
\usage{
clim.colors(n)
}
\arguments{
  \item{n}{
Number of colors to generate.
  }
}
\examples{
cols <- clim.colors(20)
lims <- seq(-1, 1, length.out = 21)
ColorBar(lims, cols)
}
\author{
History:\cr
0.0  -  2016-01  (N. Manubens, \email{nicolau.manubens at bsc.es})  -  Original code.\cr
}
\keyword{datagen}
