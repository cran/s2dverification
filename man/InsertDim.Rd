\name{InsertDim}
\alias{InsertDim}
\title{
  Adds A Dimension To A Matrix
}
\description{
  Add one dimension to the matrix var in position posdim with length lendim 
  and which correspond to lendim repetitions of the var matrix.
}
\usage{
  InsertDim(var, posdim, lendim)
}
\arguments{
  \item{var}{
    Matrix to which a dimension should be added.
  }
  \item{posdim}{
    Position of the new dimension.
  }
  \item{lendim}{
    Length of the new dimension.
  }
}
\value{
  Matrix with the added dimension.
}
\examples{
  a <- array(rnorm(15), dim = c(3, 1, 5, 1))
  print(dim(a))
  print(dim(a[, , , ]))
  print(dim(InsertDim(InsertDim(a[, , , ], 2, 1), 4, 1)))
}
\author{
  History:
  0.1  -  2011-03  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
