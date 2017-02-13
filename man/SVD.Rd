\name{SVD}
\alias{SVD}
\title{
Single Value Decomposition (Maximum Covariance Analysis)
}
\description{
Computes a Maximum Covariance Analysis (MCA) between vary and varx, both
of dimensions c(n. of time steps, n. of latitudes, n. of longitudes), each
over a region of interest, e.g.: prlr over Europe and tos over North Atlantic. 
The input fields are latitude-weighted by default (can be adjustable via 
\code{weight}).\cr
Returns a vector of squared covariance fraction (SCFs) explained by 
each pair of covariability modes, a vector of correlation coefficient 
(RUVs) between expansion coefficients (ECs) that measures their linear 
relationship, and a set of regression (MCAs) associated with the 
covariability modes (ECs). Note that MCAs are 'homogeneous' patterns obtained 
as regression/correlation between each field (predictor, predictand) 
and its expansion coefficient.\cr
The MCA is computed by default with the covariance matrix. It can be computed
with the correlation matrix by setting \code{corr = TRUE}.
}
\usage{
SVD(vary, varx, laty = NULL, latx = NULL, nmodes = 15, corr = FALSE, 
    weight = TRUE)
}
\arguments{
  \item{vary}{
Array containing the anomalies field for the predictor. The expected
dimensions are c(n. of time steps, n. of latitudes, n. of longitudes).
  }
  \item{varx}{
Array containing the anomalies field for the predictand. The expected
dimensions are c(n. of time steps, n. of latitudes, n. of longitudes).
  }
  \item{laty}{
Vector of latitudes of the array \code{vary}. Only required if 
\code{weight = TRUE}.
  }
  \item{latx}{
Vector of latitudes of the array \code{varx}. Only required if 
\code{weight = TRUE}.
  }
  \item{nmodes}{
Number of ECs/MCAs/modes retained and provided in the outputs. 
  }
  \item{corr}{
Whether to compute the MCA over a covariance matrix (FALSE) or a correlation
matrix (TRUE).
  }
  \item{weight}{
Whether to apply latitude weights on the input fields or not. TRUE by
default.
  }
}
\value{
  \item{$SC}{
Vector of squared covariance (n. of modes).
  }
  \item{$SCFs}{
Vector of squared covariance fractions (n. of modes).
  }
  \item{$RUVs}{
Vector of correlations between expansion coefficients (n. of modes).
  }
  \item{$ECs_U}{
Array of expansion coefficients of predictor field (n. of time steps, n. of modes).
  }
  \item{$MCAs_U}{
Array of covariability patterns of predictor field (c(dim), n. of modes).
  }
  \item{$ECs_V}{
Array of expansion coefficients of predictand field (n. of time steps, n. of modes).
  }
  \item{$MCAs_V}{
Array of covariability patterns of predictand field (c(dim), n. of modes).
  }
}
\examples{
# See examples on Load() to understand the first lines in this example
  \dontrun{
data_path <- system.file('sample_data', package = 's2dverification')
expA <- list(name = 'experiment', path = file.path(data_path, 
             'model/$EXP_NAME$/$STORE_FREQ$_mean/$VAR_NAME$_3hourly',
             '$VAR_NAME$_$START_DATE$.nc'))
obsX <- list(name = 'observation', path = file.path(data_path,
             '$OBS_NAME$/$STORE_FREQ$_mean/$VAR_NAME$',
             '$VAR_NAME$_$YEAR$$MONTH$.nc'))

# Now we are ready to use Load().
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- Load('tos', list(expA), list(obsX), startDates, 
                   leadtimemin = 1, leadtimemax = 4, output = 'lonlat',
                   latmin = 27, latmax = 48, lonmin = -12, lonmax = 40)
  }
  \dontshow{
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- s2dverification:::.LoadSampleData('tos', c('experiment'),
                                                c('observation'), startDates,
                                                leadtimemin = 1,
                                                leadtimemax = 4,
                                                output = 'lonlat',
                                                latmin = 27, latmax = 48,
                                                lonmin = -12, lonmax = 40)
  }
# This example computes the ECs and MCAs along forecast horizons and plots the 
# one that explains the greatest amount of variability. The example data is 
# very low resolution so it does not make a lot of sense.
ano <- Ano_CrossValid(sampleData$mod, sampleData$obs)
mca <- SVD(Mean1Dim(ano$ano_exp, 2)[1, , 1, , ], 
           Mean1Dim(ano$ano_obs, 2)[1, , 1, , ], 
           sampleData$lat, sampleData$lat)
PlotEquiMap(mca$MCAs_U[1, , ], sampleData$lon, sampleData$lat)
plot(mca$ECs_U[1, ])
PlotEquiMap(mca$MCAs_V[1, , ], sampleData$lon, sampleData$lat)
plot(mca$ECs_V[1, ])
}
\author{
History:\cr
0.1  -  2010-09  (J.-G. Serrano, \email{javier.garcia at bsc.es})  -  Original code\cr
1.0  -  2016-04  (N. Manubens, \email{nicolau.manubens at bsc.es})  -  Formatting to R CRAN
}
\keyword{datagen}
