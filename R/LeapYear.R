#'Checks Whether A Year Is Leap Year
#'
#'This function tells whether a year is a leap year or not.
#'
#'@param year A numeric value indicating the year in the Gregorian calendar.
#'
#'@return Boolean telling whether the year is a leap year or not.
#'
#'@keywords datagen
#'@author History:\cr
#'0.1  -  2011-03  (V. Guemas)  -  Original code\cr
#'1.0  -  2013-09  (N. Manubens)  -  Formatting to CRAN
#'@examples
#'print(LeapYear(1990))
#'print(LeapYear(1991))
#'print(LeapYear(1992))
#'print(LeapYear(1993))

#'@export
LeapYear <- function(year) {
  leap <- FALSE
  if (year %% 4 == 0) {
    leap <- TRUE
    if (year %% 100 == 0) {
      leap <- FALSE
      if (year %% 400 == 0) {
        leap <- TRUE
      }
    } 
  }
  
  #
  # Output
  #
  leap
}
