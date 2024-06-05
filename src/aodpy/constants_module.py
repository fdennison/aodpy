import math
from datetime import datetime, timedelta
from typing import NamedTuple

# Mathematical constants
ZERO = 0.0
ONE = 1.0
TWO = 2.0
THREE = 3.0
FOUR = 4.0
FIVE = 5.0
SIX = 6.0
SEVEN = 7.0
EIGHT = 8.0
NINE = 9.0
TEN = 10.0
HUNDRED = 100.0
THOUSAND = 1000.0
HALF = 0.5
ONE_HALF = 0.5
THREE_HALVES = 1.5
ONE_QUARTER = 0.25
THREE_QUARTERS = 0.75
ONE_THIRD = 1.0 / 3.0
TWO_THIRDS = 2.0 / 3.0

PI = 3.1415926535897932385
TWO_PI = 6.2831853071795864769
HALF_PI = 1.5707963267948966192
RADIANS_PER_DEGREE = 0.0174532925199432958
DEGREES_PER_RADIAN = 57.2957795130823208768
EULER_CONSTANT = 0.5772156649015328606

# Universal physical constants
GRAVITATIONAL_CONSTANT = 6.6705e-11  # Nm^2kg^{-2}
UNIVERSAL_GAS_CONSTANT = 8.3143435  # JK^{-1}mol^{-1}
SPEED_LIGHT = 2.9979251e+08  # m/s
BOLTZMANN_CONSTANT = 1.380546e-23  # J/K
STEFAN_BOLTZMANN_CONSTANT = 5.669710e-08  # Wm^{-2}K^{-4}
AVOGADRO_CONSTANT = 6.022529e+23  # mol^{-1}
ELECTRONIC_CHARGE = 1.602102e-19  # C
ELECTRONIC_MASS = 9.1090813e-31  # kg
PLANCK_CONSTANT = 6.625916e-34  # Js
RADIATION_CONSTANT_FIRST = 3.741509e-16  # Wm^2
RADIATION_CONSTANT_SECOND = 1.438796e-02  # mK

# Physical constants for the earth
GEO_MASS = 5.9733e+24  # kg
GEO_SEM_MAJ_AX = 6378.388e+03  # m
GEO_SEM_MIN_AX = 6356.912e+03  # m
GEO_ECC = 8.199179e-02
GEO_FLATTENING = 3.367e-03
GEO_GRAV = 1.9964985e+07  # m^{1.5}/s
GRAVITY = 9.80665  # ms^{-2}, standard value

# Physical constants for the atmosphere
SPECIFIC_HEAT_P_DRY_AIR = 1.004e+03  # J kg^{-1} K^{-1}
SPECIFIC_HEAT_V_DRY_AIR = 0.717e+03  # J kg^{-1} K^{-1}
SPECIFIC_HEAT_P_WATER_VAPOUR = 1.850e+03  # J kg^{-1} K^{-1}
SPECIFIC_HEAT_V_WATER_VAPOUR = 1.390e+03  # J kg^{-1} K^{-1}
LATENT_HEAT_VAPORIZATION = 2.500e+06  # J kg^{-1}

TRIPLE_POINT_WATER_K = 2.731600e+02
TRIPLE_POINT_WATER_C = 7.200000e-03
FREEZING_POINT_WATER_K = 2.731528e+02

# Molecular weights
MOLECULAR_WEIGHT_DRY_AIR = 2.8964e-02  # kg/mol
MOLECULAR_WEIGHT_H2O = 1.8000e-02  # kg/mol
MOLECULAR_WEIGHT_CO2 = 4.4000e-02  # kg/mol
MOLECULAR_WEIGHT_O2 = 3.2000e-02  # kg/mol
MOLECULAR_WEIGHT_O3 = 4.8000e-02  # kg/mol

# Atmospheric concentrations (default values)
VOL_MIX_RATIO_N2 = 7.8083e-01  # ppv
VOL_MIX_RATIO_O2 = 2.0947e-01  # ppv
VOL_MIX_RATIO_CO2 = 3.7000e-04  # ppv

# Time conversion factors
MINUTES_PER_DAY = 1440.0
SECONDS_PER_HOUR = 3600.0
SECONDS_PER_DAY = 86400.0
DAYS_PER_YEAR = 365.25
MILLI_SECONDS_PER_HOUR = 3600.0e+03
MILLI_SECONDS_PER_DAY = 86400.0e+03

# Length conversion factors
MICRO_METRES_PER_METRE = 1.0e+06
MICRO_METRES_PER_CENTI_METRE = 1.0e+04
NANO_METRES_PER_CENTI_METRE = 1.0e+07
NANO_METRES_PER_MICRO_METRE = 1.0e+03
SQUARE_CENTI_METRES_PER_SQUARE_METRE = 1.0e+04

# Concentration conversion factors
PARTS_PER_MILLION = 1.0e+06
PARTS_PER_BILLION = 1.0e+09
PARTS_PER_TRILLION = 1.0e+12

## time_module.py
#from datetime import datetime, timedelta
#
#class Date(NamedTuple):
#    year: int
#    month: int
#    day: int
#
#class Time(NamedTuple):
#    hour: int
#    minute: int
#    second: int
#
#class DateTime(NamedTuple):
#    year: int
#    month: int
#    day: int
#    hour: int
#    minute: int
#    second: int
#
#def days_in_month(month: int, year: int) -> int:
#    if month == 2:
#        if (year % 4 == 0 and year % 100 != 0) or year % 400 == 0:
#            return 29
#        else:
#            return 28
#    else:
#        return [31, None, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][month - 1]
#
#def day_of_year(day: int, month: int, year: int) -> int:
#    day_of_year = sum(days_in_month(m, year) for m in range(1, month))
#    return day_of_year + day
#
#def dayofyear(date: Date) -> int:
#    return day_of_year(date.day, date.month, date.year)
#
#def daybreak(year: int, numday: int) -> tuple[int, int]:
#    month = 1
#    day = 1
#    while numday > days_in_month(month, year):
#        numday -= days_in_month(month, year)
#        month += 1
#    return numday, month
#
#def num_days_year(year: int) -> int:
#    return sum(days_in_month(month, year) for month in range(1,12))
