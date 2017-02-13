\name{clim.palette}
\alias{clim.colors}
\alias{clim.palette}
\title{Generate Climate Color Palettes}
\description{
Generates a colorblind friendly color palette with color ranges useful in climate temperature variable plotting.\cr
}
\usage{
clim.palette(palette = 'bluered')

clim.colors(n, palette = 'bluered')
}
\arguments{
  \item{palette}{
Which type of palette to generate: from blue through white to red ('bluered'), from red through white to blue ('redblue'), from yellow through orange to red ('yellowred'), or from red through orange to red ('redyellow').
  }
  \item{n}{
Number of colors to generate.
  }
}
\examples{
lims <- seq(-1, 1, length.out = 21)

ColorBar(lims, color_fun = clim.palette('redyellow'))

cols <- clim.colors(20)
ColorBar(lims, cols)
}
\author{
History:\cr
0.0  -  2016-01  (N. Manubens, \email{nicolau.manubens at bsc.es})  -  Original code.\cr
}
\keyword{datagen}
