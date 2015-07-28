\name{Smoothing}
\alias{Smoothing}
\title{
Smoothes A Matrix Along A Dimension
}
\description{
Smoothes a matrix of any number of dimensions up to 10 dimensions along one of its dimensions
}
\usage{
Smoothing(var, runmeanlen = 12, numdimt = 4)
}
\arguments{
  \item{var}{
<<<<<<< HEAD
Array to be smoothed along one of its dimension (typically the forecast time dimension).
=======
Matrix to be smoothed along one of its dimension (typically the forecast time dimension).
>>>>>>> master
  }
  \item{runmeanlen}{
Running mean length in number of sampling units (typically months).
  }
  \item{numdimt}{
Dimension to smooth.
  }
}
\value{
Array with same the dimensions as 'var' but smoothed along the 'numdimt'-th dimension.
}
\examples{
# Load sample data as in Load() example:
example(Load)
clim <- Clim(sampleData$mod, sampleData$obs)
ano_exp <- Ano(sampleData$mod, clim$clim_exp)
ano_obs <- Ano(sampleData$obs, clim$clim_obs)
runmean_months <- 12
dim_to_smooth <- 4  # Smooth along lead-times
smooth_ano_exp <- Smoothing(ano_exp, runmean_months, dim_to_smooth)
smooth_ano_obs <- Smoothing(ano_obs, runmean_months, dim_to_smooth)
PlotAno(smooth_ano_exp, smooth_ano_obs, startDates, 
        toptitle = "Smoothed Mediterranean mean SST", ytitle = "K",
        fileout = "tos_smoothed_ano.eps")
}
\author{
History:\cr
0.1  -  2011-03  (V. Guemas, \email{virginie.guemas at ic3.cat})  -  Original code\cr
1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Formatting to R CRAN
}
\keyword{datagen}
