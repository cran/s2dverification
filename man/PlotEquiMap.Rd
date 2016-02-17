\name{PlotEquiMap}
\alias{PlotEquiMap}
\title{
Maps A Two-Dimensional Variable On A Cylindrical Equidistant Projection
}
\description{
Map a two dimensional matrix with (longitude, latitude) dimensions on a cylindrical equidistant latitude and longitude projection.
}
\usage{
PlotEquiMap(var, lon, lat, toptitle = "", sizetit = 1, units = "", 
            brks = NULL, cols = NULL, square = TRUE, 
            filled.continents = TRUE, contours = NULL, brks2 = NULL, 
            dots = NULL, axelab = TRUE, labW = FALSE, intylat = 20, 
            intxlon = 20, drawleg = TRUE, boxlim = NULL, 
            boxcol = 'purple2', boxlwd= 10, subsampleg = 1, 
            numbfig = 1, colNA = "white")
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
Field coloured with squares (TRUE) for each grid cell or spatial smoothing (FALSE).\cr
Default: TRUE.
  }
  \item{filled.continents}{
Continents filled in grey (TRUE) or represented by a black line (FALSE).\cr
Default = TRUE.\cr
Filling unavailable if crossing Greenwich.\cr
Filling unavailable if square = FALSE.
  }
  \item{contours}{
Matrix to be added to the plot and shown with contours. Default = NULL.
  }
  \item{brks2}{
Contour levels, optional.
  }
  \item{dots}{
Matrix with TRUE / FALSE flags to add black dots over the maps (to show where a score is significant for example). Option only available if square = TRUE.
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
Draw a colorbar. Can be FALSE only if square = FALSE.\cr
Must be FALSE if numbfig > 1. Default = TRUE.
  }
  \item{subsampleg}{
Supsampling factor of the interval between ticks on colorbar.\cr
Default = 1 = every colour level.
  }
  \item{numbfig}{
Number of figures in the final multipanel.
  }
  \item{colNA}{
Color used to represent NA. Default = 'white'
  }
  \item{boxlim}{
The area over which to draw a rectangle over the map.\cr
A vector of length 4 such that: c(western boundary, southern boundary, eastern boundary, northern boundary)
  }
  \item{boxcol}{
Color of the rectangle. Default: 'purple2'.
  }
  \item{boxlwd}{
Thickness of the rectangle side. Default: 10.
  }
}
\examples{
# See examples on Load() to understand the first lines in this example
  \dontrun{
configfile <- paste0(tempdir(), '/sample.conf')
ConfigFileCreate(configfile, confirm = FALSE)
c <- ConfigFileOpen(configfile)
c <- ConfigEditDefinition(c, 'DEFAULT_VAR_MIN', '-1e19', confirm = FALSE)
c <- ConfigEditDefinition(c, 'DEFAULT_VAR_MAX', '1e19', confirm = FALSE)
data_path <- system.file('sample_data', package = 's2dverification')
exp_data_path <- paste0(data_path, '/model/$EXP_NAME$/')
obs_data_path <- paste0(data_path, '/$OBS_NAME$/')
c <- ConfigAddEntry(c, 'experiments', dataset_name = 'experiment', 
     var_name = 'tos', main_path = exp_data_path,
     file_path = '$STORE_FREQ$_mean/$VAR_NAME$_3hourly/$VAR_NAME$_$START_DATE$.nc')
c <- ConfigAddEntry(c, 'observations', dataset_name = 'observation', 
     var_name = 'tos', main_path = obs_data_path,
     file_path = '$STORE_FREQ$_mean/$VAR_NAME$/$VAR_NAME$_$YEAR$$MONTH$.nc')
ConfigFileSave(c, configfile, confirm = FALSE)

# Now we are ready to use Load().
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- Load('tos', c('experiment'), c('observation'), startDates, 
                   output = 'lonlat', latmin = 27, latmax = 48, lonmin = -12,
                   lonmax = 40, configfile = configfile)
  }
  \dontshow{
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'), 
                                                c('observation'), startDates,
                                                output = 'lonlat', 
                                                latmin = 27, latmax = 48, 
                                                lonmin = -12, lonmax = 40) 
  }
PlotEquiMap(sampleData$mod[1, 1, 1, 1, , ], sampleData$lon, sampleData$lat, 
            toptitle = 'Predicted sea surface temperature for Nov 1960 from 1st Nov',
            sizetit = 0.5)
}
\author{
History:\cr
0.1  -  2011-11  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
0.2  -  2013-04  (R. Saurral \email{ramiro.saurral at ic3.cat}) - LabW\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to R CRAN
}
\keyword{dynamic}
