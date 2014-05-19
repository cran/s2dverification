\name{Plot2VarsVsLTime}
\alias{Plot2VarsVsLTime}
\title{
  Plot Two Scores With Confidence Intervals In A Common Plot
}
\description{
  Plots two input variables having the same dimensions in a common plot.
  One plot for all experiments.
  Input variables should have dimensions (nexp/nmod, nltime).
}
\usage{
  Plot2VarsVsLTime(var1, var2, toptitle = "", ytitle = "", monini = 1, 
                   freq = 12, nticks = NULL, limits = NULL, 
                   listexp = c("exp1", "exp2", "exp3"), 
                   listvars = c("var1", "var2"), biglab = FALSE, hlines = NULL, 
                   leg = TRUE, siglev = FALSE, sizetit = 1, 
                   fileout = "output_plot2varsvsltime.eps", show_conf = TRUE)
}
\arguments{
  \item{var1}{
    Matrix of dimensions (nexp/nmod, nltime).
  }
  \item{var2}{
    Matrix of dimensions (nexp/nmod, nltime).
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
    List of experiment names, up to three, optional.
  }
  \item{listvars}{
    List of names of input variables, optional.
  }
  \item{biglab}{
    TRUE/FALSE for presentation/paper plot. Default = FALSE.
  }
  \item{hlines}{
    c(a, b, ...) Add horizontal black lines at Y-positions a, b, ...
    Default: NULL.
  }
  \item{leg}{
    TRUE/FALSE if legend should be added or not to the plot. Default = TRUE.
  }
  \item{siglev}{
    TRUE/FALSE if significance level should replace confidence interval.
    Default = FALSE.
  }
  \item{sizetit}{
    Multiplicative factor to change title size, optional.
  }
  \item{fileout}{
    Name of output ps file.
  }
  \item{show_conf}{
    TRUE/FALSE to show/not confidence intervals for input variables.
  }
}
\details{
  Examples of input:
  ------------------
  
  RMSE error for a number of experiments and along lead-time: (nexp, nltime)
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
  required_complete_row <- 3  # Discard start dates that contain NA along lead-times
  leadtimes_per_startdate <- 60
  rms <- RMS(Mean1Dim(smooth_ano_exp, dim_to_mean), 
             Mean1Dim(smooth_ano_obs, dim_to_mean), 
             compROW = required_complete_row, 
             limits = c(ceiling((runmean_months + 1) / 2), 
                        leadtimes_per_startdate - floor(runmean_months / 2)))
  smooth_ano_exp_m_sub <- smooth_ano_exp - InsertDim(Mean1Dim(smooth_ano_exp, 2, 
                          narm = TRUE), 2, dim(smooth_ano_exp)[2])
  spread <- Spread(smooth_ano_exp_m_sub, c(2, 3))  
  Plot2VarsVsLTime(InsertDim(rms[, , , ], 1, 1), spread$sd, 
                   toptitle = 'RMSE and spread', monini = 11, freq = 12, 
                   listexp = c('CMIP5 IC3'), listvar = c('RMSE', 'spread'),
                   fileout = 'plot2vars.eps')
}
\author{
  History:
  1.0  -  2013-03  (I. Andreu-Burillo, \email{isabel.andreu-burillo@ic3.cat}) 
                                                             -  Original code
}
\keyword{dynamic}
