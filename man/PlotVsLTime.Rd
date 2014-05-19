\name{PlotVsLTime}
\alias{PlotVsLTime}
\title{
  Plots A Score Along The Forecast Time With Its Confidence Interval
}
\description{
  Plots The Correlation (\code{Corr()}) or the Root Mean Square Error 
  (\code{RMS()}) between the forecasted values and their observational 
  counterpart or the slopes of their trends (\code{Trend()}) or the 
  InterQuartile Range, Maximum-Mininum, Standard Deviation or Median 
  Absolute Deviation of the Ensemble Members (\code{Spread()}), or
  the ratio between the Ensemble Spread and the RMSE of the Ensemble 
  Mean (\code{RatioSDRMS()}) along the forecast time for all the input 
  experiments on the same figure with their confidence intervals.
}
\usage{
  PlotVsLTime(var, toptitle = "", ytitle = "", monini = 1, freq = 12, 
  nticks = NULL, limits = NULL, listexp = c("exp1", "exp2", "exp3"), 
  listobs = c("obs1", "obs2", "obs3"), biglab = FALSE, hlines = NULL, leg = TRUE, 
  siglev = FALSE, fileout = "output_plotvsltime.eps", sizetit = 1, show_conf = TRUE)
}
\arguments{
  \item{var}{
    Matrix containing any Prediction Score with dimensions:
         (nexp/nmod, 3/4 ,nltime)
      or (nexp/nmod, nobs, 3/4 ,nltime)
  }
  \item{toptitle}{
    Main title, optional.
  }
  \item{ytitle}{
    Title of Y-axis, optional.
  }
  \item{monini}{
    Starting month between 1 and 12. Default = 1.
  }
  \item{freq}{
    1 = yearly, 12 = monthly, 4 = seasonal, ... Default = 12.
  }
  \item{nticks}{
    Number of ticks and labels on the x-axis, optional.
  }
  \item{limits}{
    c(lower limit, upper limit): limits of the Y-axis, optional.
  }
  \item{listexp}{
    List of experiment names, optional.
  }
  \item{listobs}{
    List of observation names, optional.
  }
  \item{biglab}{
    TRUE/FALSE for presentation/paper plot. Default = FALSE.
  }
  \item{hlines}{
    c(a,b, ..) Add horizontal black lines at Y-positions a,b, ...
    Default = NULL.
  }
  \item{leg}{
    TRUE/FALSE if legend should be added or not to the plot. Default = TRUE.
  }
  \item{siglev}{
    TRUE/FALSE if significance level should replace confidence interval.
    Default = FALSE.
  }
  \item{fileout}{
    Name of the output eps file.
  }
  \item{sizetit}{
    Multiplicative factor to change title size, optional.
  }
  \item{show_conf}{
    TRUE/FALSE to show/not confidence intervals for input variables.
  }
}
\details{
  Examples of input:


     Model and observed output from \code{Load()} then \code{Clim()} then
     \code{Ano()} then \code{Smoothing()}:
       (nmod, nmemb, nsdate, nltime) and (nobs, nmemb, nsdate, nltime)
     then averaged over the members 
       \code{Mean1Dim(var_exp/var_obs, posdim = 2)}:
       (nmod, nsdate, nltime) and (nobs, nsdate, nltime)
     then passed through 
       \code{Corr(exp, obs, posloop = 1, poscor = 2)} or 
       \code{RMS(exp, obs, posloop = 1, posRMS = 2)}:
       (nmod, nobs, 3, nltime)
     would plot the correlations or RMS between each exp & each obs as a
     function of the forecast time.
}
\examples{
  startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
  sampleData <- Load('tos', c('i00k'), c('ERSST'), startDates, nleadtime = 124, 
                     leadtimemin = 1, leadtimemax = 60, output = 'areave', 
                     latmin = 30, latmax = 45, lonmin = 0, lonmax = 40)
  clim <- Clim(sampleData$mod, sampleData$obs)
  ano_exp <- Ano(sampleData$mod, clim$clim_exp)
  ano_obs <- Ano(sampleData$obs, clim$clim_obs)
  runmean_months <- 12
  dim_to_smooth <- 4  # Smooth along lead-times
  smooth_ano_exp <- Smoothing(ano_exp, runmean_months, dim_to_smooth)
  smooth_ano_obs <- Smoothing(ano_obs, runmean_months, dim_to_smooth)
  dim_to_mean <- 2  # Mean along members
  required_complete_row <- 3  # Discard startdates for which there are NA leadtimes
  leadtimes_per_startdate <- 60
  corr <- Corr(Mean1Dim(smooth_ano_exp, dim_to_mean), 
               Mean1Dim(smooth_ano_obs, dim_to_mean), 
               compROW = required_complete_row, 
               limits = c(ceiling((runmean_months + 1) / 2), 
                          leadtimes_per_startdate - floor(runmean_months / 2)))
  PlotVsLTime(corr, toptitle = "correlations", ytitle = "correlation", 
              monini = 11, limits = c(-1, 2), listexp = c('CMIP5 IC3'), 
              listobs = c('ERSST'), biglab = FALSE, hlines = c(-1, 0, 1), 
              fileout = 'tos_cor.eps')
}
\author{
  History:
  0.1  -  2011-03  (V. Guemas,                           -  Original code
                 \email{virginie.guemas@ic3.cat})
  0.2  -  2013-03  (I. Andreu-Burillo,                   -  Introduced parameter
       -         \email{isabel.andreu-burillo@ic3.cat})  -  sizetit
  0.3  -  2013-10  (I. Andreu-Burillo,                   -  Introduced parameter
       -         \email{isabel.andreu-burillo@ic3.cat})  -  show_conf 
  1.0  -  2013-11  (N. Manubens,
                 \email{nicolau.manubens@ic3.cat})       -  Formatting to CRAN
}
\keyword{dynamic}
