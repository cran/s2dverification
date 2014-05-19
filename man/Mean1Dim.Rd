\name{Mean1Dim}
\alias{Mean1Dim}
\title{
  Averages A Matrix Along A Dimension
}
\description{
  Averages the matrix var along the posdim dimension between limits [1] and 
  limits [2] if limits argument is provided by the user.
}
\usage{
  Mean1Dim(var, posdim, narm = TRUE, limits = NULL)
}
\arguments{
  \item{var}{
    Matrix to average.
  }
  \item{posdim}{
    Dimension to average along.
  }
  \item{narm}{
    Ignore NA (TRUE) values or not (FALSE).
  }
  \item{limits}{
    Limits to average between.
  }
}
\value{
  Matrix with one dimension less than the input one containing the average along posdim dimension.
}
\examples{
  a <- array(rnorm(24), dim = c(2, 3, 4))
  print(a)
  print(Mean1Dim(a, 2))
}
\author{
  History:
  0.1  -  2011-04  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
