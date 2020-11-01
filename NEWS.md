# s2dverification 2.9.0 (Release date: 2020-10-30)
- Correct the time retrieval in Load() when start date and the first lead time in netCDF file do not match. 
In addition, when the lead time in each data set is not consistent, the retrieved data should change according to
 the first data set in the same Load call.
- One new parameter 'K' is added in Composite() to indicate the maximum number of composites.
- Revise the per-paired method in Clim() when NA exists.
- Correct the Corr() p-value.
- Bugfix for CDO version reading. The error occurred when the CDO version number is followed by letters.
- Bugfix for Ano() when obs and exp have inconsistent dimensions.

# s2dverification 2.8.6 (Release date: 2019-10-17)
- Apply Roxygen2 format to all the files.
- Bug fix in Composite().
- Bug fix in Ano(). Recommend to assign the dimensions by name to avoid confusion when the dimensions have same length.
- Trend() documentation error fix.
- Introduce new function PlotMatrix().

# s2dverification 2.8.0 (Release date: 2017-02-13)
- Change licence from GPL-3 to LGPL-3.
- New veriApply compatible score functions (.BrierScore, .Corr, .RMS, .RMSSS, .RatioRMS, .RatioSDRMS and .Trend).
- New function CDORemap to interpolate R arrays with CDO.
- New function ArrayToNetCDF to save R arrays with metadata in NetCDF files.
- Enhance plot_timeseries.R and plot_maps.R example scripts to work with file-per-chunk data, for auto-ecearth v3.2.1a.
- Add colour-blind colour bars for the map plots.
- Add warning in Load when extrapolating data.
- Bug fix in ProbBins when called with cross-validation.
- Enhance documentation.
- Adapt UltimateBrier to SpecsVerification 0.5-0.
- Enhancements to adjust size and resolution in plotting functions.
- Solve PlotEquiMap bug when values equal to the lower limit.
- Bug fix in Ano.
- Bug fix in PlotVsLTime.
- Small update in the configuration file.

# s2dverification 2.7.0 (Release date: 2016-08-24)
Enhanced PlotEquiMap() and PlotStereoMap() with lots of new options and fixed issues:
- Colour bar with triangle ends and lots of new features.
- Margins, labels, ticks, colour bar, titles now arranged properly.
- Now possibile to adjust colour and border of continents, size, colour and type of contour lines, size of labels, ticks and margins, colour and width of boxes, etc.
- Draw multiple superimposed dot/symbol layers.
- Draw boxes in PlotStereoMap().
- PlotStereoMap() with bounding circle.
- Add function PlotLayout() to automatically generate complex figure layouts and populate with plots from multi-dimensional arrays.
- Fix and updated corrupted example scripts (required for new auto-ecearth releases to work).
- Add function Subset() to easily take chunks of data arrays.
- Fix bug in Load() under some particular configurations.
- Enhance margins in PlotAno().
- Update sample data to be together with metadata as provided by Load().
- Update and fix in the BSC Load() configuration file. 

# s2dverification 2.6.0 (Release date: 2016-06-06)
- Update configuration file.
- Functions to compute variability modes and project data on EOF() and ProjectField().
- Function to compute co·variability modes: SVD(), by Javi.
- Function to compute the NAO: NAO(), by Fabian, Virginie, Lauriane, Martin.
- Brier score/skill score accounting for small ensemble/start date sample size: UltimateBrier().
- K-means spatial clustering: Cluster().
- Synthetic data generator: ToyModel().
- Tropical cyclone downscaling: StatSeasAtlHurr().
- Function to composite fields: Composite().
- Function to generate map animations: AnimateMap().
- Function to plot time series with box-and-whisker plots: PlotBoxWhisker().
- Possible to disable computation of confidence intervals or p-values in ACC(), Corr(), RatioRMS(), RatioSDRMS(), RMS(), RMSSS(), Spread() and Trend().
- Possible to adjust confidence level in all functions that provide confidence intervals: ACC(), Corr(), RMS(), Spread() and Trend().
- Possible to plot arrows in PlotEquiMap().
- Possible to save plots in multiple formats, to file or onscreen from all plot functions.
- Objects returned by Load() are now closer to the format in downscaleR. The initial and end date of each time step is provided now.
- Enhancements in Smoothing().
- Load() now stops if the tag $START_DATE$/$YEAR$+$MONTH$+$DAY$ is not in the path pattern of any of the experiments.

# s2dverification 2.5.0 (Release date: 2016-01-04)
- Fix bugs when using masks in Load()
- Able to specify masks with paths to NetCDF files

# s2dverification 2.4.7 (Release date: 2015-11-15)
- Update plot_timeseries.R to new paths and to 'ncdf4'.
- Improve performance when retrieving subsets of data (regions of earth or time periods).
- Add possibility to use Load() without a configuration file. See details on parameters 'exp' and 'obs' in ?Load  

Load() now returns plenty of metadata. Highlighted:  
- Paths to all loaded files
- Paths to not found files
- Stamp with all the provided parameters to exactly reproduce any Load() call
- The name of the common grid (if any), following CDO's naming conventions  

Other enhancements in Load():
- Enhance error handling and error messages
- Add “progress bar”
- Detect automatically grid of the files. No need to specify it
- Detect automatically if the requested variable is 2-dimensional or global mean. No need to specify it
- Possibility to load observations only, from a limited period of time only
- Possibility to load NetCDF files with disordered dimensions
- Remove system dependency of 's2dverification' to NCO and some GNU tools
- Simplify configuration file: removed lists of variables and reduced from 5 tables to 2, one for experimental datasets and one for observational datasets. You can convert old configuration files to the new format with the script in /shared/earth/software/scripts/convertConfig.R as follows: /shared/earth/software/scripts/convertConfig.R /path/to/configfile.conf
- Fix and updated the sample script plot_timeseries.R
- Fix wrong entries in BSC configuration file for some ice variables.

# s2dverification 2.4.0 (Release date: 2015-07-27)
- Option to draw rectangles in PlotEquiMap()
- Motification of Corr() to accomodate ranked correlation
- Add the possibility to load the second set of HadCM3 decadal data (i3p1)
- Add functions to assist in manipulating the configuration file
- Improve examples that use extremely reduced experimental and observational datasets
- Uniformize documentation style
- Add possibility to configure dimension names to look for inside NetCDF files
- Add the possibility to load ESA observations from SMHI
- Fix bug that happened in some cases when a common grid is not specified

# s2dverification 2.3.2 (Release date: 2015-04-23)
- New CRPS() function to compute the continuous ranked probability score for ensemble forecasts.
- New ProbBins() function to compute probabilistic information of a forecast relative to a threshold or a quantile.
- Load() stops and warns if the masks provided are not in the correct grid.
- Load() didn't apply, as expected, the same masks in observations as in experiments when possible. Now fixed.
- Enhancement in Clim() documentation.
- Enhancements in Load() and configuration file documentation.
- HadSLP dataset is now loadable

# s2dverification 2.3.1 (Release date: 2015-03-09)
- Loading observations only fixed
- Loading only one leadtime fixed
- Loading a 2D variable when the first observation was not stored in file-per-dataset format fixed
- Parameter 'ncores' changed to 'nprocs'
- Improvements in configuration file mechanism and documentation

# s2dverification 2.3.0 (Release date: 2015-03-02)
- Configuration file mechanism to specify new dataset or variable paths, grids, etc.
- New parameters in Load() to specify maximum and minimum values.
- New supported dataset formats. See '?Load' in R after loadings2dverificationfor more information.
- More efficient memory usage in Load() and usage of multiple parallel processes (faster).
- NetCDF4 + OPeNDAP support

# s2dverification 2.2.0 (Release date: 2014-12-16)
- ACC provides confidence intervals obtained with bootstrap method
- Function to plot ACC score
- Function to plot variables on a polar stereographic projection
- Possibility of loading observations only
- Possibility to load more ice variables
- Adjustable significance level in the Corr function
- Adjustable number size in ColorBar

# s2dverification 2.1.0 (Release date: 2014-01-23)
- Demo scripts 'plot_timeseries.R' and 'plot_maps.R' available in the 'inst/doc' directory in thes2dverificationrepository.
- Documentation on how to specify the grids and masks to the function Load() has been added to its help page, code and package manual.

# s2dverification 2.0.0 (Release date: 2013-08-02)
- Use of the standard R package structure.
- Use of the google's R style guide.
- Functions that involved RClim set of funcions have been kept apartfrom the package (AnimVsLTime, BlueRed, PlotMap, ProjMap) as well as the authors.
- New functions have been added: Alpha, EnoNew, Filter, FitAcfCoef, FitAutocor, GenSeries, Spectrum.
- Extended help.
