\name{ConfigShowTable}
\alias{ConfigShowTable}
\alias{ConfigShowDefinitions}
\title{
Show Configuration Tables And Definitions
}
\description{
These functions show the tables of supported datasets and definitions in a configuration object obtained via ConfigFileCreate() or ConfigFileOpen().
}
\usage{
ConfigShowTable(configuration, dataset_type, line_numbers = NULL)
ConfigShowDefinitions(configuration)
}
\arguments{
  \item{configuration}{
Configuration object obtained from ConfigFileCreate() or ConfigFileOpen().
  }
  \item{dataset_type}{
In ConfigShowTable(), 'dataset_type' tells whether the table to show is of experimental datasets or of observational datasets.\cr
Can take values 'experiments' or 'observations'.
  }
  \item{line_numbers}{
'line_numbers' is an optional vector of numbers as long as the number of entries in the specified table.\cr
Intended for internal use.
  }
}
\value{
These functions return nothing.
}
\seealso{ConfigApplyMatchingEntries, ConfigEditDefinition, ConfigEditEntry, ConfigFileOpen, ConfigShowSimilarEntries, ConfigShowTable}
\author{
History:\cr
0.1  -  2015-05  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  First version
1.0  -  2015-11  (N. Manubens, \email{nicolau.manubens at ic3.cat})  -  Removed grid column and storage formats
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
configuration <- ConfigAddEntry(configuration, "experiments", "last", 
                 "ExampleExperiment2", "ExampleVariable", 
                 "/path/to/ExampleExperiment2/", 
                 "ExampleVariable/ExampleVariable_$START_DATE$.nc")
# Edit entry to generalize for any variable. Changing variable needs .
configuration <- ConfigEditEntry(configuration, "experiments", 1, 
                 var_name = ".*", 
                 file_path = "$VAR_NAME$/$VAR_NAME$_$START_DATE$.nc")
# Show tables, lists and definitions
ConfigShowTable(configuration, 'experiments')
ConfigShowDefinitions(configuration)
}
\keyword{datagen}
