\name{Spread}
\alias{Spread}
\title{
  Computes InterQuartile Range, Maximum-Minimum, Standard Deviation and Median
  Absolute Deviation of the Ensemble Members
}
\description{
  Computes the InterQuartile Range, the Maximum minus Mininum, the Standard 
  Deviation and the Median Absolute Deviation along the list of dimensions
  provided by the posdim argument (typically along the ensemble member 
  and start date dimension).
  The confidence interval is computed by bootstrapping.
}
\usage{
  Spread(var, posdim = 2, narm = TRUE)
}
\arguments{
  \item{var}{
    Matrix of any number of dimensions up to 10.
  }
  \item{posdim}{
    List of dimensions along which to compute IQR/MaxMin/SD/MAD.
  }
  \item{narm}{
    TRUE/FALSE if NA removed/kept for computation. Default = TRUE.
  }
}
\details{
  Example:
  --------
  To compute IQR, Max-Min, SD & MAD accross the members and start dates of var 
  output from \code{Load()} or \code{Ano()} or \code{Ano_CrossValid()}, call:
    spread(var, posdim = c(2, 3), narm = TRUE)
}
\value{
  Matrix with the same dimensions as var except along the first posdim
  dimension which is replaced by a length 3 dimension, corresponding to the 
  lower limit of the 95\% confidence interval, the spread and the upper limit 
  of the 95\% confidence interval for each experiment/leadtime/latitude/longitude.
  \item{$iqr}{
    InterQuartile Range.
  }
  \item{$maxmin}{
    Maximum - Minimum.
  }
  \item{$sd}{
    Standard Deviation.
  }
  \item{$mad}{
    Median Absolute Deviation.
  }
}
\examples{
  startDates <- c('19901101', '19951101', '20001101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'areave', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  clim <- Clim(sampleData$mod, sampleData$obs)
  ano_exp <- Ano(sampleData$mod, clim$clim_exp)
  runmean_months <- 12
  dim_to_smooth <- 4  # Smooth along lead-times
  smooth_ano_exp <- Smoothing(ano_exp, runmean_months, dim_to_smooth)
  smooth_ano_exp_m_sub <- smooth_ano_exp - InsertDim(Mean1Dim(smooth_ano_exp, 2, 
                          narm = TRUE), 2, dim(smooth_ano_exp)[2])
  spread <- Spread(smooth_ano_exp_m_sub, c(2, 3))
  PlotVsLTime(spread$iqr, 
              toptitle = "Inter-Quartile Range between ensemble members",
              ytitle = "K", monini = 11, limits = NULL, 
              listexp = c('CMIP5 IC3'), listobs = c('ERSST'), biglab = FALSE, 
              hlines = c(0), fileout = 'tos_iqr.eps')
  PlotVsLTime(spread$maxmin, toptitle = "Maximum minus minimum of the members", 
              ytitle = "K", monini = 11, limits = NULL, 
              listexp = c('CMIP5 IC3'), listobs = c('ERSST'), biglab = FALSE, 
              hlines = c(0), fileout = 'tos_maxmin.eps')
  PlotVsLTime(spread$sd, toptitle = "Standard deviation of the members", 
              ytitle = "K", monini = 11, limits = NULL, 
              listexp = c('CMIP5 IC3'), listobs = c('ERSST'), biglab = FALSE, 
              hlines = c(0), fileout = 'tos_sd.eps')
  PlotVsLTime(spread$mad, toptitle = "Median Absolute Deviation of the members",
              ytitle = "K", monini = 11, limits = NULL, 
              listexp = c('CMIP5 IC3'), listobs = c('ERSST'), biglab = FALSE, 
              hlines = c(0), fileout = 'tos_mad.eps')
}
\author{
  History:
  0.1  -  2011-03  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{datagen}
