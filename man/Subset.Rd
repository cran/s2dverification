% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Subset.R
\name{Subset}
\alias{Subset}
\title{Subset a Data Array}
\usage{
Subset(x, along, indices, drop = FALSE)
}
\arguments{
\item{x}{A multidimensional array to be sliced. It can have dimension names 
either in \code{names(dim(x))} or either in the attribute 'dimensions'.}

\item{along}{Vector with references to the dimensions to take the subset 
from: either integers or dimension names.}

\item{indices}{List of indices to take from each dimension specified in 
'along'. If a single dimension is specified in 'along' the indices can be 
directly provided as a single integer or as a vector.}

\item{drop}{Whether to drop all the dimensions of length 1 in the resulting 
array, none, only those that are specified in 'along', or only those that 
are not specified in 'along'. The possible values are, respectively: 'all' 
or TRUE, 'none' or FALSE, 'selected', and 'non-selected'.}
}
\description{
This function allows to subset (i.e. slice, take a chunk of) an array, in a 
similar way as done in the function \code{take()} in the package plyr. There
are two main inprovements:\cr\cr
The input array can have dimension names, either 
in \code{names(dim(x))} or in the attribute 'dimensions', and the dimensions 
to subset along can be specified via the parameter \code{along} either with 
integer indices or either by their name.\cr\cr
There are additional ways to adjust which dimensions are dropped in the 
resulting array: either to drop all, to drop none, to drop only the ones that 
have been sliced or to drop only the ones that have not been sliced.\cr\cr
If an array is provided without dimension names, dimension names taken from 
the parameter \code{dim_names} will be added to the array.
}
\examples{
subset <- Subset(sampleMap$mod, c('dataset', 'sdate', 'ftime'), 
                list(1, 1, 1), drop = 'selected')
PlotLayout(PlotEquiMap, c('lat', 'lon'), subset, 
          sampleMap$lon, sampleMap$lat, 
          titles = paste('Member', 1:3))

}
\keyword{datagen}
