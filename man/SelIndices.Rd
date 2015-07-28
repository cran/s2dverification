\name{SelIndices}
\alias{SelIndices}
\title{
Slices A Matrix Along A Dimension
}
\description{
This function allows to select a subensemble from a matrix of any dimensions, providing the dimension along which the user aims at cutting the input matrix and between which indices.
}
\usage{
SelIndices(var, posdim, limits)
}
\arguments{
  \item{var}{
A matrix of any dimensions.
  }
  \item{posdim}{
The dimension along which a submatrix should be selected.
  }
  \item{limits}{
The lower and upper indice of the selection along the posdim dimension.
  }
}
\value{
The sliced matrix.
}
\examples{
a <- array(rnorm(24), dim = c(2, 3, 4, 1))
print(a)
print(a[, , 2:3, ])
print(dim(a[, , 2:3, ]))
print(SelIndices(a, 3, c(2, 3)))
print(dim(SelIndices(a, 3, c(2, 3))))
}
\author{
History:\cr
0.1  -  2011-04  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
