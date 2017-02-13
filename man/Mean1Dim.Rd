\name{Mean1Dim}
\alias{Mean1Dim}
\title{
Averages An Array Along A Dimension
}
\description{
Averages the array along the posdim dimension along the user specified dimension. The user can specify a subset of the dimension to take the mean along.
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
Limits to average between. Default is to take the mean along the entire dimension.
  }
}
\value{
Array with one dimension less than the input array, containing the average along the posdim dimension.
}
\examples{
a <- array(rnorm(24), dim = c(2, 3, 4))
print(a)
print(Mean1Dim(a, 2))
}
\author{
History:\cr
0.1  -  2011-04  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
