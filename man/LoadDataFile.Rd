\name{.LoadDataFile}
\alias{.LoadDataFile}
\title{
Load Data From File Into Environment
}
\description{
This function receives a 'work piece', a named list which contains information on a data file to be loaded.\cr
It can be run in 'explore' mode or in 'load' mode.\cr
When running in 'explore' mode, the metadata of the file is read and the sizes of the dimensions in the file are returned in a named list:
  \itemize{
    \item{'nmemb': Number of members}
    \item{'nltime': Number of lead-times}
    \item{'lon': Longitudes in the file}
    \item{'lat': Latitudes in the file}
  }
When running in 'load' mode it loads and performs any requested computations in additional parameters in the 'work piece', such as interpolating, slicing, ..., and finally stores it in a shared memory matrix pointed by the parameter 'out_pointer' in the 'work piece'.
}
\usage{
.LoadDataFile(work_piece, explore_dims = FALSE, silent = FALSE)
}
\arguments{
  \item{work_piece}{
Named list with information on the file to load or explore and additional parameters.\cr
The needed variables in the work piece are:
    \itemize{
      \item{In explore mode:
        \itemize{
          \item{dataset_type: 'exp'/'obs'}
          \item{filename: full path to the data file}
          \item{dimnames: named list with names 'longitudes', 'latitudes' and 'members' and with values the associated actual dimension names inside the NetCDF files.}
          \item{namevar: name of the variable in the nc file}
          \item{is_2d_var: TRUE/FALSE}
          \item{grid: common grid into which interpolate or NULL if no interpolation}
          \item{remap: interpolation method ('remapbil'/'remapdis'/'remapcon'/'remapbic')}
          \item{lon_limits: c(lon_min, lon_max)}
          \item{lat_limits: c(lat_min, lat_max)}
          \item{is_file_per_member: TRUE/FALSE}
        }
      }
      \item{In load mode:
        \itemize{
          \item{dataset_type: 'exp'/'obs'}
          \item{filename: full path to the data file}
          \item{dimnames: named list with names 'longitudes', 'latitudes' and 'members' and with values the associated actual dimension names inside the NetCDF files.}
          \item{namevar: name of the variable in the nc file}
          \item{is_2d_var: TRUE/FALSE}
          \item{grid: common grid into which interpolate or NULL if no interpolation}
          \item{remap: interpolation method ('remapbil'/'remapdis'/'remapcon'/'remapbic')}
          \item{lon_limits: c(lon_min, lon_max)}
          \item{lat_limits: c(lat_min, lat_max)}
          \item{is_file_per_dataset: TRUE/FALSE}
          \item{startdates: in the case that filename points to a whole dataset file, the list of starting dates to load must be specified. c('sdate1', 'sdate2', ...)}
          \item{out_pointer: big.matrix descriptor pointing to the array (transformed into a matrix) where to keep the data}
          \item{dims: named list with dimension sizes of the original array into which the data is kept. Names must be c(['nlon', ]['nlat', ]'nltime', 'nmember', 'nsdates', 'ndat').}
          \item{indices: vector of initial positions corresponding to each dimension in 'dims' where to store data in the original array. First two indices ('nlon', 'nlat') can be missing.}
          \item{nmember: number of members expected to be loaded}
          \item{mask: complete (untrimmed + interpolated if needed) mask to activate/deactivate data points, with dimensions c('lon', 'lat')}
          \item{leadtimes: vector of time indices to be loaded from the file}
          \item{var_limits: c(var_min, var_max)}
        }
      }
    }
  }
  \item{explore_dims}{
Run in dimension explore mode (TRUE) or in load and calculation mode (FALSE).\cr
Takes by default the value FALSE (calculation mode).
  }
  \item{silent}{
Parameter to allow (FALSE) or deactivate (TRUE) printing of explanatory messages.\cr
When deactivated any warning messages will still be displayed.\cr
Takes by default the value FALSE (verbose mode).
  }
}
\value{
When called in 'explore' mode, a named list is returned with the found lengths for members, leadtimes and the latitudes and longitudes already trimmed and reordered if needed. The names are 'nmemb', 'nltime', 'lon', 'lat'. If it is a file from a file-per-member dataset, the number of files that match the filename replacing the $MEMBER$ part by an asterisk is returned (which is the supposed number of members). There are known issues with this method of detection. See documentation on parameter 'nmember' and 'nmemberobs' in Load() function.\cr
When called in 'calculation' mode, TRUE is returned if the file was found and FALSE is returned otherwise.  
}
\examples{
data <- s2dverification:::.LoadDataFile(list(dataset_type = 'exp', 
        filename = system.file('sample_data/model/experiment/monthly_mean', 
                               'tos_3hourly/tos_19901101.nc', 
                               package = 's2dverification'),
        namevar = 'tos', is_2d_var = TRUE, lon_limits = c(-12, 40), 
        lat_limits = c(27, 48), is_file_per_member = FALSE, 
        dimnames = list(longitudes = 'longitude', latitudes = 'latitude', 
        members = 'ensemble')), explore_dims = TRUE, silent = FALSE)
}
\author{
History:\cr
0.1  -  2015-01  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  First version
}
\keyword{datagen}
