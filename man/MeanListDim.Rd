\name{MeanListDim}
\alias{MeanListDim}
\title{
Averages An Array Along Multiple Dimensions
}
\description{
Averages an array along a set of dimensions given by the argument dims.
}
\usage{
MeanListDim(var, dims, narm = TRUE)
}
\arguments{
  \item{var}{
Input array.
  }
  \item{dims}{
List of dimensions to average along.
  }
  \item{narm}{
Ignore NA (TRUE) values or not (FALSE).
  }
}
\value{
The averaged array, with the dimensions specified in \code{dims} removed.
}
\examples{
a <- array(rnorm(24), dim = c(2, 3, 4))
print(a)
print(Mean1Dim(a, 2))
print(MeanListDim(a, c(2, 3)))
}
\author{
History:\cr
0.1  -  2011-04  (V. Guemas, \email{vguemas@ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to R CRAN\cr
1.1  -  2015-03  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Improved memory usage
}
\keyword{datagen}
