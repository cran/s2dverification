\name{MeanListDim}
\alias{MeanListDim}
\title{
Averages A Matrix Along Various Dimensions
}
\description{
Averages the matrix var along a set of dimensions given by the argument dims.
}
\usage{
MeanListDim(var, dims, narm = TRUE)
}
\arguments{
  \item{var}{
Matrix to average.
  }
  \item{dims}{
List of dimensions to average along.
  }
  \item{narm}{
Ignore NA (TRUE) values or not (FALSE).
  }
}
\value{
Matrix with as many dimensions less than the input matrix as provided by the list dims and containing the average along this list of dimensions.
}
\examples{
a <- array(rnorm(24), dim = c(2, 3, 4))
print(a)
print(Mean1Dim(a, 2))
print(MeanListDim(a, c(2, 3)))
}
\author{
History:\cr
0.1  -  2011-04  (V. Guemas, \email{vguemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
