\name{PlotEquiMap}
\alias{PlotEquiMap}
\title{
  Maps A Two-Dimensional Variable On A Cylindrical Equidistant Projection
}
\description{
  Map a two dimensional matrix with (longitude, latitude) dimensions on 
  a cylindrical equidistant latitude and longitude projection.
}
\usage{
  PlotEquiMap(var, lon, lat, toptitle = "", sizetit = 1, units = "", 
              brks = NULL, cols = NULL, square = TRUE, 
              filled.continents = TRUE, contours = NULL, brks2 = NULL, 
              dots = NULL, axelab = TRUE, labW = FALSE, intylat = 20, intxlon = 20, 
              drawleg = TRUE, subsampleg = 1, numbfig = 1, colNA = "white")
}
\arguments{
  \item{var}{
    Matrix to plot with (longitude, latitude) dimensions.
  }
  \item{lon}{
    Array of longitudes.
  }
  \item{lat}{
    Array of latitudes.
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
  \item{square}{
    Field coloured with squares (TRUE) for each grid cell or spatial smoothing 
    (FALSE). 
    Default: TRUE.
  }
  \item{filled.continents}{
    Continents filled in grey (TRUE) or represented by a black line (FALSE). 
    Default = TRUE.
    Filling unavailable if crossing Greenwich.
    Filling unavailable if square = FALSE.
  }
  \item{contours}{
    Matrix to be added to the plot and shown with contours. Default = NULL.
  }
  \item{brks2}{
    Contour levels, optional.
  }
  \item{dots}{
    Matrix with TRUE / FALSE flags to add black dots over the maps (to show 
    where a score is significant for example). Option only available if 
    square = TRUE.
  }
  \item{axelab}{
    TRUE/FALSE, label the axis. Default = TRUE.
  }
  \item{labW}{
    Label the longitude axis with W instead of minus. Default = FALSE. 
  }
  \item{intylat}{
    Interval between latitude ticks on y-axis. Default = 20deg.
  }
  \item{intxlon}{
    Interval between longitude ticks on x-axis. Default = 20deg.
  }
  \item{drawleg}{
    Draw a colorbar. Can be FALSE only if square = FALSE.
    Must be FALSE if numbfig > 1. Default = TRUE.
  }
  \item{subsampleg}{
    Supsampling factor of the interval between ticks on colorbar.
    Default = 1 = every colour level.
  }
  \item{numbfig}{
    Number of figures in the final multipanel.
  }
  \item{colNA}{
    Color used to represent NA. Default = 'white'
  }
}
\examples{
  startDates <- c('19901101', '19951101', '20001101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'lonlat', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  PlotEquiMap(sampleData$mod[1, 1, 1, 1, , ], sampleData$lon, sampleData$lat, 
              toptitle = 'Predicted sea surface temperature for Nov 1960 from 1st Nov',
              sizetit = 0.5)
}
\author{
  History:
  0.1  -  2011-11  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  0.2  -  2013-04  (R. Saurral \email{ramiro.saurral@ic3.cat}) - LabW
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to R CRAN
}
\keyword{dynamic}
