\name{PlotStereoMap}
\alias{PlotStereoMap}
\title{
Maps A Two-Dimensional Variable On A Polar Stereographic Projection
}
\description{
Map a two dimensional matrix with (longitude, latitude) dimensions on a polar stereographic projection.
}
\usage{
PlotStereoMap(var, lon, lat, latlims = c(60,90), toptitle = "", sizetit = 1, 
             units = "", brks = NULL, cols = NULL, filled.continents = FALSE, 
             dots = NULL, intlat = 10, drawleg = TRUE, subsampleg = 1, 
             colNA = "white")
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
  \item{latlims}{
Latitudinal limits of the figure.\cr
Example : c(60, 90) for the North Pole\cr
          c(-90,-60) for the South Pole
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
  \item{filled.continents}{
Continents filled in grey (TRUE) or represented by a black line (FALSE).\cr
Default = FALSE.
  }
  \item{dots}{
Matrix with TRUE / FALSE flags to add black dots over the maps (to show where a score is significant for example).  
  }
  \item{intlat}{
Interval between latitude circles. Default = 10deg.
  }
  \item{drawleg}{
Draw a colorbar. 
  }
  \item{subsampleg}{
Supsampling factor of the interval between ticks on colorbar.\cr
Default = 1 = every colour level.
  }
  \item{colNA}{
Color used to represent NA. Default = 'white'
  }
}
\author{
History:\cr
1.0  -  2014-07  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code
}
\keyword{dynamic}
