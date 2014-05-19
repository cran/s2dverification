\name{PlotACC}
\alias{PlotACC}
\title{
  Plot Plumes/Timeseries Of Anomaly Correlation Coefficients
}
\description{
  Plots plumes/timeseries of ACC from a matrix with dimensions (output from
    \code{ACC()}): 
    c(nexp, nobs, nsdates, nltime, 4)
  with the fourth dimension of length 4 containing the lower limit of the 95\% 
  confidence interval, the ACC, the upper limit of the 95\% confidence 
  interval and the 95\% significance level given by a one-sided T-test.
}
\usage{
  PlotACC(ACC, sdates, toptitle = "", sizetit = 1, ytitle = "", limits = NULL, 
          legends = NULL, freq = 12, biglab = FALSE, fill = FALSE, linezero = FALSE, 
          points = TRUE, vlines = NULL, fileout = "output_PlotACC.eps")
}
\arguments{
  \item{ACC}{
    ACC matrix with with dimensions:
      c(nexp, nobs, nsdates, nltime, 4)
    with the fourth dimension of length 4 containing the lower limit of the 95\% 
    confidence interval, the ACC, the upper limit of the 95\% 
    confidence interval and the 95\% significance level.
  }
  \item{sdates}{
    List of startdates: c('YYYYMMDD','YYYYMMDD').
  }
  \item{toptitle}{
    Main title, optional.
  }
  \item{sizetit}{
    Multiplicative factor to scale title size, optional.
  }
  \item{ytitle}{
    Title of Y-axis for each experiment: c('',''), optional.
  }
  \item{limits}{
    c(lower limit, upper limit): limits of the Y-axis, optional.
  }
  \item{legends}{
    List of flags (characters) to be written in the legend, optional.
  }
  \item{freq}{
    1 = yearly, 12 = monthly, 4 = seasonal, ... Default: 12.
  }
  \item{biglab}{
    TRUE/FALSE for presentation/paper plot, Default = FALSE.
  }
  \item{fill}{
    TRUE/FALSE if filled confidence interval. Default = FALSE.
  }
  \item{linezero}{
    TRUE/FALSE if a line at y=0 should be added. Default = FALSE.
  }
  \item{points}{
    TRUE/FALSE if points instead of lines. Default = TRUE.
  }
  \item{vlines}{
    List of x location where to add vertical black lines, optional.
  }
  \item{fileout}{
    Name of the output eps file.
  }
}
\examples{
  startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 4, output = 'lonlat', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  sampleData$mod <- Season(sampleData$mod, 4, 11, 12, 2)
  sampleData$obs <- Season(sampleData$obs, 4, 11, 12, 2)
  clim <- Clim(sampleData$mod, sampleData$obs)
  ano_exp <- Ano(sampleData$mod, clim$clim_exp)
  ano_obs <- Ano(sampleData$obs, clim$clim_obs)
  acc <- ACC(Mean1Dim(sampleData$mod, 2), 
             Mean1Dim(sampleData$obs, 2))
  PlotACC(acc$ACC, startDates, toptitle = "Anomaly Correlation Coefficient")
}
\author{
  History:
  0.1  -  2013-08  (V. Guemas, \email{virginie.guemas@ic3.cat})  -  Original code
  1.0  -  2013-09  (N. Manubens, \email{nicolau.manubens@ic3.cat})  -  Formatting to CRAN
}
\keyword{dynamic}
