% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConfigEditEntry.R
\name{ConfigEditEntry}
\alias{ConfigEditEntry}
\alias{ConfigAddEntry}
\alias{ConfigRemoveEntry}
\title{Add, Remove Or Edit Entries In The Configuration}
\usage{
ConfigEditEntry(
  configuration,
  dataset_type,
  position,
  dataset_name = NULL,
  var_name = NULL,
  main_path = NULL,
  file_path = NULL,
  nc_var_name = NULL,
  suffix = NULL,
  varmin = NULL,
  varmax = NULL
)

ConfigAddEntry(
  configuration,
  dataset_type,
  position = "last",
  dataset_name = ".*",
  var_name = ".*",
  main_path = "*",
  file_path = "*",
  nc_var_name = "*",
  suffix = "*",
  varmin = "*",
  varmax = "*"
)

ConfigRemoveEntry(
  configuration,
  dataset_type,
  dataset_name = NULL,
  var_name = NULL,
  position = NULL
)
}
\arguments{
\item{configuration}{Configuration object obtained via ConfigFileOpen() 
or ConfigFileCreate() that will be modified accordingly.}

\item{dataset_type}{Whether to modify a table of experimental datasets or 
a table of observational datasets. Can take values 'experiments' or 
'observations' respectively.}

\item{position}{'position' tells the index in the table of the entry to 
edit or remove. Use ConfigShowTable() to see the index of the entry.
In ConfigAddEntry() it can also take the value "last" (default), that will 
put the entry at the end of the corresponding level, or "first" at the 
beginning. See ?ConfigFileOpen for more information.
If 'dataset_name' and 'var_name' are specified this argument is ignored in 
ConfigRemoveEntry().}

\item{dataset_name, var_name, main_path, file_path, nc_var_name, suffix, varmin, varmax}{These parameters tell the dataset name, variable name, main path, ..., of 
the entry to add, edit or remove.\cr 'dataset_name' and 'var_name' can take 
as a value a POSIX 1003.2 regular expression (see ?ConfigFileOpen).\cr
Other parameters can take as a value a shell globbing expression 
(see ?ConfigFileOpen).\cr
'dataset_name' and 'var_name' take by default the regular expression '.*' 
(match any dataset and variable name), and the others take by default '*' 
(associate to the pair 'dataset_name' and 'var_name' all the defined 
default values. In this case '*' has a special behaviour, it won't be 
used as a shell globbing expression. See ?ConfigFileOpen and 
?ConfigShowDefinitions).\cr
'var_min' and 'var_max' must be a character string.\cr
To define these values, you can use defined variables via $VARIABLE_NAME$ 
or other entry attributes via $ATTRIBUTE_NAME$. See ?ConfigFileOpen for 
more information.}
}
\value{
The function returns an accordingly modified configuration object. 
 To apply the changes in the configuration file it must be saved using 
 ConfigFileSave().
}
\description{
ConfigAddEntry(), ConfigEditEntry() and ConfigRemoveEntry() are functions 
to manage entries in a configuration object created with ConfigFileOpen().\cr
Before adding an entry, make sure the defaults don't do already what you 
want (ConfigShowDefinitions(), ConfigShowTable()).\cr
Before adding an entry, make sure it doesn't override and spoil what other 
entries do (ConfigShowTable(), ConfigFileOpen()).\cr
Before adding an entry, make sure there aren't other entries that already 
do what you want (ConfigShowSimilarEntries()).
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
                "last", "ExampleExperiment", "ExampleVariable", 
                "/path/to/ExampleExperiment/", 
                "ExampleVariable/ExampleVariable_$START_DATE$.nc")
# Add another entry
configuration <- ConfigAddEntry(configuration, "experiments",
                "last", "ExampleExperiment2", "ExampleVariable", 
                "/path/to/ExampleExperiment2/",
                "ExampleVariable/ExampleVariable_$START_DATE$.nc")
# Edit second entry to generalize for any variable. Changing variable needs .
configuration <- ConfigEditEntry(configuration, "experiments", 2, 
                var_name = ".*",
                file_path = "$VAR_NAME$/$VAR_NAME$_$START_DATE$.nc")
# Remove first entry
configuration <- ConfigRemoveEntry(configuration, "experiments",
                "ExampleExperiment", "ExampleVariable")
# Show results
ConfigShowTable(configuration, "experiments")
# Save the configuration
ConfigFileSave(configuration, config_file, confirm = FALSE)
}
\seealso{
ConfigApplyMatchingEntries, ConfigEditDefinition, ConfigEditEntry, 
 ConfigFileOpen, ConfigShowSimilarEntries, ConfigShowTable
}
\author{
History:\cr
 0.1 - 2015-05 (N. Manubens) - First version\cr
 1.0  - 2015-11 (N. Manubens) - Removed grid column and storage formats
}
\keyword{datagen}
