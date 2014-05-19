\name{sampleMap}
\docType{data}
\alias{sampleMap}
\title{Sample of Observational and Experimental Data for Forecast Verification
       In Function Of Longitudes And Latitudes}
\description{
  This data set provides data in function of longitudes and latitudes for the 
  variable 'tos', i.e. sea surface temperature, from the decadal climate 
  prediction experiment run at IC3 in the context of the CMIP5 project. 
  Its name within IC3 local database is 'i00k'. 
  The observational dataset used for verification is the 'ERSST' observational 
  dataset, in this example.
  
  The data is provided through a variable named 'sampleMap' and is structured 
  as expected from the 'Load()' function in the 's2dverification' package if 
  was called as follows:
  
    sampleData <- Load('tos', c('i00k'), c('ERSST'),  c('19851101', '19901101', 
                       '19951101', '20001101', '20051101'), nleadtime = 124, 
                       leadtimemin = 1, leadtimemax = 60, sampleperiod = 1, 
                       output = 'lonlat', latmin = 30, latmax = 45, lonmin = 0, 
                       lonmax = 40) 

  Check the documentation on 'Load()' in the package 's2dverification' for 
  more information.
}
\usage{
  data(sampleMap)
}
\format{
  The data set provides with a variable named 'sampleMap'.

  sampleMap$mod is an array that contains the experimental data and the 
  dimension meanings and values are:
    c(# of experimental datasets, # of members, # of starting dates, 
      # of lead-times, # of latitudes, # of longitudes)
    c(1, 5, 5, 60, 13, 36)

  sampleMap$obs is an array that contains the observational data and the
  dimension meanings and values are:
    c(# of observational datasets, # of members, # of starting dates,
      # of lead-times, # of latitudes, # of longitudes)
    c(1, 1, 5, 60, 13, 36)

  sampleMap$lat is an array with the 13 latitudes covered by the data.

  sampleMap$lon is an array with the 36 longitudes covered by the data.
}
\keyword{datasets}
