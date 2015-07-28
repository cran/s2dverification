\name{ConfigEditEntry}
\alias{ConfigAddEntry}
\alias{ConfigEditEntry}
\alias{ConfigRemoveEntry}
\title{
Add, Remove Or Edit Entries In The Configuration
}
\description{
ConfigAddEntry(), ConfigEditEntry() and ConfigRemoveEntry() are functions to manage entries in a configuration object created with ConfigFileOpen().\cr
Before adding an entry, make sure the defaults don't do already what you want (ConfigShowDefinitions(), ConfigShowTable()).\cr
Before adding an entry, make sure it doesn't override and spoil what other entries do (ConfigShowTable(), ConfigFileOpen()).\cr
Before adding an entry, make sure there aren't other entries that already do what you want (ConfigShowSimilarEntries()).\cr
}
\usage{
ConfigAddEntry(configuration, dataset_type, store_format, position = "last", 
               dataset_name = ".*", var_name = ".*", main_path = "*", 
               file_path = "*", nc_var_name = "*", grid = "*", suffix = "*", 
               varmin = "*", varmax = "*")
ConfigEditEntry(configuration, dataset_type, store_format, position, 
                dataset_name = NULL, var_name = NULL, main_path = NULL, 
                file_path = NULL, nc_var_name = NULL, grid = NULL, 
                suffix = NULL, varmin = NULL, varmax = NULL)
ConfigRemoveEntry(configuration, dataset_type, store_format, 
                  dataset_name = NULL, var_name = NULL, position = NULL)
}
\arguments{
  \item{configuration}{
Configuration object obtained via ConfigFileOpen() or ConfigFileCreate() that will be modified accordingly.
  }
  \item{dataset_type}{
Whether to modify a table of experimental datasets or a table of observational datasets. Can take values 'experiments' or 'observations' respectively.
  }
  \item{store_format}{
There's a table associated to each storage format (see ?ConfigFileOpen).\cr
'store_format' tells which one to load. Can take values 'file-per-startdate', 'file-per-member' (if 'dataset_type' is 'experiments'), or 'file-per-member', 'file-per-month', 'file-per-dataset' (if 'dataset_type' is 'observations').
  }
  \item{position}{
'position' tells the index in the table of the entry to edit or remove. Use ConfigShowTable() to see the index of the entry.\cr
In ConfigAddEntry() it can also take the value "last" (default), that will put the entry at the end of the corresponding level, or "first" at the beginning. See ?ConfigFileOpen for more information.\cr
If 'dataset_name' and 'var_name' are specified this argument is ignored in ConfigRemoveEntry().
  }
  \item{dataset_name,var_name,main_path,file_path,nc_var_name,grid,suffix,varmin,varmax}{
These parameters tell the dataset name, variable name, main path, ..., of the entry to add, edit or remove.\cr
'dataset_name' and 'var_name' can take as a value a POSIX 1003.2 regular expression (see ?ConfigFileOpen).\cr
Other parameters can take as a value a shell globbing expression (see ?ConfigFileOpen).\cr
'dataset_name' and 'var_name' take by default the regular expression '.*' (match any dataset and variable name), and the others take by default '*' (associate to the pair 'dataset_name' and 'var_name' all the defined default values. In this case '*' has a special behaviour, it won't be used as a shell globbing expression. See ?ConfigFileOpen and ?ConfigShowDefinitions).\cr
'var_min' and 'var_max' must be a character string.\cr
To define these values, you can use defined variables via $VARIABLE_NAME$ or other entry attributes via $ATTRIBUTE_NAME$. See ?ConfigFileOpen for more information.\cr
  }
}
\value{
The function returns an accordingly modified configuration object. To apply the changes in the configuration file it must be saved using ConfigFileSave().
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
# table which will associate the experiment "ExampleExperiment" and variable 
# "ExampleVariable" to some information about its location.
configuration <- ConfigAddEntry(configuration, "experiments", 
                 "file-per-startdate", "last", "ExampleExperiment", 
                 "ExampleVariable", "/path/to/ExampleExperiment", 
                 "ExampleVariable/ExampleVariable_$START_DATE$.nc")
# Add another entry
configuration <- ConfigAddEntry(configuration, "experiments",
                 "file-per-startdate", "last", "ExampleExperiment2",
                 "ExampleVariable", "/path/to/ExampleExperiment2",
                 "ExampleVariable/ExampleVariable_$START_DATE$.nc")
# Edit second entry to generalize for any variable. Changing variable needs .
configuration <- ConfigEditEntry(configuration, "experiments",
                 "file-per-startdate", 2, var_name = ".*",
                 file_path = "$VAR_NAME$/$VAR_NAME$_$START_DATE$.nc")
# Remove first entry
configuration <- ConfigRemoveEntry(configuration, "experiments",
                 "file-per-startdate", "ExampleExperiment", "ExampleVariable")
# Show results
ConfigShowTable(configuration, "experiments", "file-per-ensemble")
# Now put variable in list of supported variables. In this example will be a 
# two-dimensional variable.
configuration <- ConfigAddVar(configuration, "2d", "ExampleMistake")
# Look in list of variables and see there's nothing similar to what we look 
# for. The previous introduction was mistaken
configuration <- ConfigShowSimilarVars(configuration, "Variable")
# Remove variable and add the desired one
configuration <- ConfigRemoveVar(configuration, "2d", "ExampleMistake")
configuration <- ConfigAddVar(configuration, "2d", "ExampleVariable")
# See results
ConfigShowVars(configuration)
# Save the configuration
ConfigFileSave(configuration, config_file, confirm = FALSE)
}
\keyword{datagen}
