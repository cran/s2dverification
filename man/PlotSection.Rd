\name{PlotSection}
\alias{PlotSection}
\title{
  Plots A Vertical Section
}
\description{
  Plot a (longitude,depth) or (latitude,depth) section.
}
\usage{
  PlotSection(var, horiz, depth, toptitle = "", sizetit = 1, units = "", 
              brks = NULL, cols = NULL, axelab = TRUE, intydep = 200, 
              intxhoriz = 20, drawleg = TRUE)
}
\arguments{
  \item{var}{
    Matrix to plot with (longitude/latitude, depth) dimensions.
  }
  \item{horiz}{
    Array of longitudes or latitudes.
  }
  \item{depth}{
    Array of depths.
  }
  \item{toptitle}{
    Title, optional.
  }
  \item{sizetit}{
    Multiplicative factor to increase title size, optional.
  }
  \item{units}{
    Units, optional.
  }
  \item{brks}{
    Colour levels, optional.
  }
  \item{cols}{
    List of colours, optional.
  }
  \item{axelab}{
    TRUE/FALSE, label the axis. Default = TRUE.
  }
  \item{intydep}{
    Interval between depth ticks on y-axis. Default: 200m.
  }
  \item{intxhoriz}{
    Interval between longitude/latitude ticks on x-axis. 
    Default: 20deg.
  }
  \item{drawleg}{
    Draw colorbar. Default: TRUE.
  }
}
\examples{
  sampleData <- sampleDepthData
  PlotSection(sampleData$mod[1, 1, 1, 1, , ], sampleData$lat, sampleData$depth,
              toptitle = 'temperature 1995-11 member 0')
}
\author{
  History:
  0.1  -  2012-09  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{dynamic}
