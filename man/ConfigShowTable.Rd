\name{ConfigShowTable}
\alias{ConfigShowTable}
\alias{ConfigShowVars}
\alias{ConfigShowDefinitions}
\title{
Show Configuration Tables Lists And Definitions
}
\description{
These functions show the tables of supported datasets, lists of supported variables and definitions in a configuration object obtained via ConfigFileCreate() or ConfigFileOpen().
}
\usage{
ConfigShowTable(configuration, dataset_type, store_format, line_numbers = NULL)
ConfigShowDefinitions(configuration)
ConfigShowVars(configuration)
}
\arguments{
  \item{configuration}{
Configuration object obtained from ConfigFileCreate() or ConfigFileOpen().
  }
  \item{dataset_type}{
In ConfigShowTable(), 'dataset_type' tells whether the table to show is of experimental datasets or of observational datasets.\cr
Can take values 'experiments' or 'observations'.
  }
  \item{store_format}{
In ConfigShowTable(), 'store_format' tells the store format of the table to show.\cr
Can take values 'file-per-startdate', 'file-per-member' if 'dataset_type' is 'experiments'.\cr
Can take values 'file-per-month', 'file-per-member', 'file-per-dataset' if 'dataset_type' is 'observations'.\cr
  }
  \item{line_numbers}{
'line_numbers' is an optional vector of numbers as long as the number of entries in the specified table.\cr
Intended for internal use.
  }
}
\value{
These functions return nothing.
}
\seealso{ConfigAddVar, ConfigApplyMatchingEntries, ConfigEditDefinition, ConfigEditEntry, ConfigFileOpen, ConfigMatchingVars, ConfigShowSimilarEntries, ConfigShowTable}
\author{
History:\cr
0.1  -  2015-05  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  First version
}
\examples{
# Create an empty configuration file
config_file <- paste0(tempdir(), "/example.conf")
ConfigFileCreate(config_file, confirm = FALSE)
# Open it into a configuration object
configuration <- ConfigFileOpen(config_file)
# Add an entry at the bottom of 4th level of file-per-startdate experiments 
# table which will associate the experiment "ExampleExperiment2" and variable 
# "ExampleVariable" to some information about its location.
configuration <- ConfigAddEntry(configuration, "experiments", 
                 "file-per-startdate", "last", "ExampleExperiment2", 
                 "ExampleVariable", "/path/to/ExampleExperiment2", 
                 "ExampleVariable/ExampleVariable_$START_DATE$.nc")
# Edit entry to generalize for any variable. Changing variable needs .
configuration <- ConfigEditEntry(configuration, "experiments", 
                 "file-per-startdate", 1, var_name = ".*", 
                 file_path = "$VAR_NAME$/$VAR_NAME$_$START_DATE$.nc")
# Now put a pattern to the list of supported variables. 
# In this example will be a two-dimensional variable.
configuration <- ConfigAddVar(configuration, "2d", "Example.*Variable")
# Show tables, lists and definitions
ConfigShowTable(configuration, 'experiments', 'file-per-startdate')
ConfigShowVars(configuration)
ConfigShowDefinitions(configuration)
}
\keyword{datagen}
