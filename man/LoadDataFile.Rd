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
    \item{'member': Number of members}
    \item{'time': Number of lead-times}
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
          \item{dimnames: named list with names 'lon', 'lat' and 'member' and with values the associated actual dimension names inside the NetCDF files.}
          \item{namevar: name of the variable in the nc file}
          \item{grid: common grid into which interpolate or NULL if no interpolation}
          \item{remap: interpolation method ('remapbil'/'remapdis'/'remapcon'/'remapbic')}
          \item{lon_limits: c(lon_min, lon_max)}
          \item{lat_limits: c(lat_min, lat_max)}
          \item{is_file_per_member: TRUE/FALSE}
          \item{is_file_per_dataset: TRUE/FALSE}
          \item{is_single_dataset: TRUE/FALSE}
        }
      }
      \item{In load mode:
        \itemize{
          \item{dataset_type: 'exp'/'obs'}
          \item{filename: full path to the data file}
          \item{dimnames: named list with names 'lon', 'lat' and 'member' and with values the associated actual dimension names inside the NetCDF files.}
          \item{namevar: name of the variable in the nc file}
          \item{is_2d_var: TRUE/FALSE}
          \item{grid: common grid into which interpolate or NULL if no interpolation}
          \item{remap: interpolation method ('remapbil'/'remapdis'/'remapcon'/'remapbic')}
          \item{lon_limits: c(lon_min, lon_max)}
          \item{lat_limits: c(lat_min, lat_max)}
          \item{is_file_per_dataset: TRUE/FALSE}
          \item{startdates: in the case that filename points to a whole dataset file, the list of starting dates to load must be specified. c('sdate1', 'sdate2', ...)}
          \item{out_pointer: big.matrix descriptor pointing to the array (transformed into a matrix) where to keep the data}
          \item{dims: named list with dimension sizes of the original array into which the data is kept. Names must be c(['lon', ]['lat', ]'time', 'member', 'sdates', 'dataset').}
          \item{indices: vector of initial positions corresponding to each dimension in 'dims' where to store data in the original array. First two indices ('lon', 'lat') can be missing.}
          \item{nmember: number of members expected to be loaded}
          \item{mask: complete (untrimmed + interpolated if needed) mask to activate/deactivate data points, with dimensions c('lon', 'lat'), see \code{?Load} for more details.}
          \item{leadtimes: vector of time indices to be loaded from the file}
          \item{var_limits: c(var_min, var_max)}
          \item{is_single_dataset: TRUE/FALSE (whether the user only asked for data of one single dataset. If so, then data won't be re-interpolated when the first longitude of its grid is != 0. Otherwise it must be re-interpolated to ensure all data will be properly aligned over longitudes.}
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
  \item{remapcells}{
Width in number of cells of the surrounding area of the requested subset to be taken into account for the interpolation. See parameter \code{remapcells} of \code{Load()}.
  }
}
\value{
When called in 'explore' mode, a named list is returned with:
  \itemize{
    \item{dims}{
List with the found lengths for members, leadtimes and the latitudes and longitudes already trimmed and reordered if needed. The names are 'member', 'time', 'lon', 'lat'. If it is a file from a file-per-member dataset, the number of files that match the filename replacing the $MEMBER_NUMBER$ part by an asterisk is returned (which is the supposed number of members). There are known issues with this method of detection. See documentation on parameter 'nmember' and 'nmemberobs' in Load() function. When the specified file path is a URL, the returned number of members is NULL.\cr
    }
    \item{is_2d_var}{
Boolean indicating whether the found variable is 2-dimensional (TRUE) or a global mean (FALSE).
    }
    \item{grid}{
Character string with the name of the grid of the file, following the cdo grid naming conventions.
    }
    \item{var_long_name}{
A character string with the variable long name. If not available, the short name is returned.
    }
    \item{units}{
A character string with the units of the variable.
    }
  }
When called in 'calculation' mode, the found file path or URL is returned if the file was found and NULL is returned otherwise.  
}
\examples{
  \dontrun{
data <- s2dverification:::.LoadDataFile(list(dataset_type = 'exp', 
        filename = system.file('sample_data/model/experiment/monthly_mean', 
                               'tos_3hourly/tos_19901101.nc', 
                               package = 's2dverification'),
        namevar = 'tos', lon_limits = c(-12, 40), 
        lat_limits = c(27, 48), is_file_per_member = TRUE, 
        dimnames = list(lon = 'longitude', lat = 'latitude', 
        member = 'ensemble')), explore_dims = TRUE, silent = FALSE)
  }
}
\author{
History:\cr
0.1  -  2015-01  (N. Manubens)  -  First version
}
\keyword{datagen}
\keyword{internal}
