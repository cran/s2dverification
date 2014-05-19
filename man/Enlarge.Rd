\name{Enlarge}
\alias{Enlarge}
\title{
  Extends The Number Of Dimensions of A Matrix
}
\description{
  Extends the number of dimensions of var to numdims (the added dimensions have 
  length 1).
}
\usage{
  Enlarge(var, numdims)
}
\arguments{
  \item{var}{
    Matrix to be extended.
  }
  \item{numdims}{
    Output number of dimensions.
  }
}
\value{
  Extended matrix.
}
\examples{
  data <- array(1, c(2, 2, 3))
  print(dim(Enlarge(data, 5)))
}
\author{
  History:
  0.1  -  2011-03  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
