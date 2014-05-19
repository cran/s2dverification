\name{Regression}
\alias{Regression}
\title{
  Computes The Regression Of A Matrix On Another Along A Dimension
}
\description{
  Computes the regression of the input matrice vary on the input matrice varx 
  along the posREG dimension by least square fitting. Provides the slope of 
  the regression, the associated confidence interval, and the intercept.
  Provides also the vary data filtered out from the regression onto varx.
  The confidence interval relies on a student-T distribution.
}
\usage{
  Regression(vary, varx, posREG = 2)
}
\arguments{
  \item{vary}{
    Matrix of any number of dimensions up to 10.
  }
  \item{varx}{
    Matrix of any number of dimensions up to 10. Same dimensions as vary.
  }
  \item{posREG}{
    Position along which to compute the regression.
  }
}
\value{
  \item{$regression}{
    Matrix with same dimensions as varx and vary except along posREG dimension 
    which is replaced by a length 4 dimension, corresponding to the lower limit 
    of the 95\% confidence interval, the slope, the upper limit of the 95\% 
    confidence interval and the intercept.
  }
  \item{$filtered}{
     Same dimensions as vary filtered out from the regression onto varx along 
     the posREG dimension.
  }
}
\examples{
  startDates <- c('19901101', '19951101', '20001101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 4, output = 'lonlat', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  sampleData$mod <- Season(sampleData$mod, 4, 11, 12, 2)
  sampleData$obs <- Season(sampleData$obs, 4, 11, 12, 2)
  reg <- Regression(Mean1Dim(sampleData$mod, 2),
                    Mean1Dim(sampleData$obs, 2), 2)
  PlotEquiMap(reg$regression[1, 2, 1, , ], sampleData$lon, sampleData$lat, 
              toptitle='Regression of the prediction on the observations', 
              sizetit = 0.5)
}
\author{
  History:
  0.1  -  2013-05  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
