\name{Eno}
\alias{Eno}
\title{
  Computes Effective Sample Size With Classical Method
}
\description{
  Computes the effective number of independant data along the posdim 
  dimension of a matrix. This effective number of independant date may be 
  required to perform statistical/inference tests.
  Based on eno function from Caio Coelho from rclim.txt.
}
\usage{
  Eno(obs, posdim)
}
\arguments{
  \item{obs}{
    Matrix of any number of dimensions up to 10.
  }
  \item{posdim}{
    Dimension along which to compute the effective sample size.
  }
}
\value{
  Same dimensions as var but without the posdim dimension.
}
\examples{
  startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'lonlat', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  sampleData$mod <- Season(sampleData$mod, 4, 11, 1, 12)
  eno <- Eno(sampleData$mod[1, 1, , 1, , ], 1)
  PlotEquiMap(eno, sampleData$lon, sampleData$lat)
}
\author{
  History:
  0.1  -  2011-05  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
