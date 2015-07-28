\name{ConfigEditDefinition}
\alias{ConfigEditDefinition}
\alias{ConfigRemoveDefinition}
\title{
Add Modify Or Remove Variable Definitions In Configuration
}
\description{
These functions help in adding, modifying or removing variable definitions in a configuration object obtained wit ConfigFileOpen() or ConfigFileCreate().\cr
ConfigEditDefinition() will add the definition if not existing.
}
\usage{
ConfigEditDefinition(configuration, name, value, confirm = TRUE)
ConfigRemoveDefinition(configuration, name)
}
\arguments{
  \item{configuration}{
Configuration object obtained wit ConfigFileOpen() or ConfigFileCreate().
  }
  \item{name}{
Name of the variable to add/modify/remove.
  }
  \item{value}{
Value to associate to the variable.
  }
  \item{confirm}{
Flag to stipulate whether to ask for confirmation if the variable is being modified. Takes by default TRUE.
  }
}
\value{
A modified configuration object is returned.
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
# Now put variable in list of supported variables. In this example will be a 
# two-dimensional variable.
configuration <- ConfigAddVar(configuration, "2d", "ExampleVariable")
# Change the default grid
configuration <- ConfigEditDefinition(configuration, "DEFAULT_GRID", 
                                      "r256x128", confirm = FALSE)
# Now apply matching entries for variable and experiment name and show the 
# result
match_info <- ConfigApplyMatchingEntries(configuration, 'tas', 
              exp = c('ExampleExperiment2'), show_result = TRUE)
}
\keyword{datagen}
