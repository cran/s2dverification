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
            intxhoriz = 20, drawleg = TRUE, fileout = NULL, 
            width = 8, height = 5, size_units = 'in', res = 100, ...)
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
Interval between longitude/latitude ticks on x-axis.\cr
Default: 20deg.
  }
  \item{drawleg}{
Draw colorbar. Default: TRUE.
  }
  \item{width}{
File width, in the units specified in the parameter size_units (inches by default). Takes 8 by default.
  }
  \item{height}{
File height, in the units specified in the parameter size_units (inches by default). Takes 5 by default.
  }
  \item{size_units}{
Units of the size of the device (file or window) to plot in. Inches ('in') by default. See ?Devices and the creator function of the corresponding device.
  }
  \item{res}{ 
Resolution of the device (file or window) to plot in. See ?Devices and the creator function of the corresponding device.
  }
  \item{fileout}{
Name of output file. Extensions allowed: eps/ps, jpeg, png, pdf, bmp 
and tiff. \cr
Default = NULL
  }
  \item{...}{
    Arguments to be passed to the method. Only accepts the following
    graphical parameters:

    adj ann ask bg bty cex.lab cex.sub cin col.axis col.lab col.main col.sub cra crt csi cxy err family fg fig fin font font.axis font.lab font.main font.sub lend lheight ljoin lmitre lty lwd mex mfcol mfrow mfg mkh oma omd omi page pch pin plt pty smo srt tcl usr xaxp xaxs xaxt xlog xpd yaxp yaxs yaxt ylbias ylog

    For more information about the parameters see `par`
  }
}
\examples{
sampleData <- s2dverification::sampleDepthData
PlotSection(sampleData$mod[1, 1, 1, 1, , ], sampleData$lat, sampleData$depth,
            toptitle = 'temperature 1995-11 member 0')
}
\author{
History:\cr
0.1  -  2012-09  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{dynamic}
