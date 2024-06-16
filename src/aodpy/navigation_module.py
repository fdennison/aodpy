import math
import datetime as dt
from math import pi, sin, cos, tan, atan, atan2, sqrt, radians, degrees

# Define constants
LONG_REAL = float
PI = math.pi
TWO_PI = 2 * math.pi
GEO_SEMI_MAJ_AX = 6378.137  # Earth's semi-major axis (km)
GEO_SEMI_MIN_AX = 6356.7523  # Earth's semi-minor axis (km)
GEO_ECC = (GEO_SEMI_MAJ_AX - GEO_SEMI_MIN_AX) / GEO_SEMI_MAJ_AX  # Earth's eccentricity
RADIANS_PER_DEGREE = PI / 180.0
SECONDS_PER_DAY = 86400.0
GEO_RADIUS = 6378137.0
GEO_ECC = 0.08181919
GEO_SEMIMAJAX = 6378137.0
J_2 = 1082.28e-6
GEO_GRAV = 3.986005e14

# Define data types
class OrbitalElements:
    def __init__(self, id_code=0, semi_maj_ax=7.0e6, ecc=0.0, mean_anom=0.0,
                 inc=0.0, arg_perigee=0.0, asc_node=0.0, period=0.0, freq=0.0,
                 anom_period=0.0, anom_freq=0.0, prec_perigee=0.0, prec_asc_node=0.0,
                 epoch=(2000, 1, 0, 0, 0, 0, 0), jul=0, day=0.0):
        self.id_code = id_code
        self.semi_maj_ax = semi_maj_ax
        self.ecc = ecc
        self.mean_anom = mean_anom
        self.inc = inc
        self.arg_perigee = arg_perigee
        self.asc_node = asc_node
        self.period = period
        self.freq = freq
        self.anom_period = anom_period
        self.anom_freq = anom_freq
        self.prec_perigee = prec_perigee
        self.prec_asc_node = prec_asc_node
        self.epoch = epoch
        self.jul = jul
        self.day = day

class SatelliteAxes:
    def __init__(self, r=0.0, s=0.0, t=(0.0, 0.0, 0.0), u=(0.0, 0.0, 0.0),
                 v=(0.0, 0.0, 0.0), w=(0.0, 0.0, 0.0)):
        self.r = r
        self.s = s
        self.t = t
        self.u = u
        self.v = v
        self.w = w

#class OrbitalElements(NamedTuple):
#    mean_anom: float
#    inc: float
#    arg_perigee: float
#    asc_node: float
#    semi_maj_ax: float
#    ecc: float
#    anom_period: float
#    anom_freq: float
#    prec_perigee: float
#    prec_asc_node: float
#    jul: int
#    day: float
#    epoch: float

# Define default orbital elements
default_elements = OrbitalElements()
greenwich_elements = OrbitalElements(
    semi_maj_ax=1.0, period=0.997269566340, anom_period=0.997269566340,
    #epoch=(1988, 1, 0, 17, 21, 35, 354)
    epoch=dt.datetime(1988,1,1,17,21,35,354) #-dt.timedelta(days=1)
)
sun_elements = OrbitalElements(
    id_code=1, semi_maj_ax=1.496e11, ecc=0.167133900e-01, mean_anom=0.621622191e+01,
    inc=0.409120047e+00, arg_perigee=0.493460338e+01, period=0.365259635e+03,
    anom_period=0.365259635e+03, prec_perigee=0.000000822e+00, 
    #epoch=(1988, 1, 0, 0, 0, 0, 0)
    epoch=dt.datetime(1988,1,1,0,0,0,0) #-dt.timedelta(days=1)
)

GEO_FACTOR = (GEO_ECC / (1 - GEO_ECC)) * (GEO_ECC / (1 + GEO_ECC))

def geo_radius(geo_pos):
    return GEO_SEMI_MAJ_AX / math.sqrt(1 + GEO_FACTOR * geo_pos[2] ** 2)

def hit_planet(rad_sat, geo_sat, sat_fov):
    missed = False
    len_sat_fov = 0.0
    len_geo_fov = 0.0
    geo_fov = (0.0, 0.0, 0.0)

    s = (GEO_SEMI_MAJ_AX / GEO_SEMI_MIN_AX) ** 2
    t = (GEO_SEMI_MAJ_AX / rad_sat) ** 2
    p = geo_sat[0] ** 2 + geo_sat[1] ** 2 + geo_sat[2] ** 2 * s
    q = sat_fov[0] ** 2 + sat_fov[1] ** 2 + sat_fov[2] ** 2 * s
    r = geo_sat[0] * sat_fov[0] + geo_sat[1] * sat_fov[1] + geo_sat[2] * sat_fov[2] * s
    discriminant = r ** 2 - q * (p - t)

    if r >= 0 or discriminant <= 0:
        missed = True
    else:
        missed = False
        len_sat_fov = rad_sat * (abs(r) - math.sqrt(discriminant)) / q
        geo_fov = (rad_sat * geo_sat[0] + len_sat_fov * sat_fov[0],
                   rad_sat * geo_sat[1] + len_sat_fov * sat_fov[1],
                   rad_sat * geo_sat[2] + len_sat_fov * sat_fov[2])
        len_geo_fov = math.sqrt(geo_fov[0] ** 2 + geo_fov[1] ** 2 + geo_fov[2] ** 2)
        geo_fov = (geo_fov[0] / len_geo_fov, geo_fov[1] / len_geo_fov, geo_fov[2] / len_geo_fov)

    return missed, len_sat_fov, len_geo_fov, geo_fov

def hit_planet_tunnel(rad_sat, geo_sat, sat_fov):
    missed = False
    len_sat_fov = 0.0
    len_geo_fov = 0.0
    geo_fov = (0.0, 0.0, 0.0)

    s = (GEO_SEMI_MAJ_AX / GEO_SEMI_MIN_AX) ** 2
    t = (GEO_SEMI_MAJ_AX / rad_sat) ** 2
    p = geo_sat[0] ** 2 + geo_sat[1] ** 2 + geo_sat[2] ** 2 * s
    q = sat_fov[0] ** 2 + sat_fov[1] ** 2 + sat_fov[2] ** 2 * s
    r = geo_sat[0] * sat_fov[0] + geo_sat[1] * sat_fov[1] + geo_sat[2] * sat_fov[2] * s
    discriminant = r ** 2 - q * (p - t)

    if r >= 0 or discriminant <= 0:
        missed = True
    else:
        missed = False
        len_sat_fov = rad_sat * (abs(r) + math.sqrt(discriminant)) / q
        geo_fov = (rad_sat * geo_sat[0] + len_sat_fov * sat_fov[0],
                   rad_sat * geo_sat[1] + len_sat_fov * sat_fov[1],
                   rad_sat * geo_sat[2] + len_sat_fov * sat_fov[2])
        len_geo_fov = math.sqrt(geo_fov[0] ** 2 + geo_fov[1] ** 2 + geo_fov[2] ** 2)
        geo_fov = (geo_fov[0] / len_geo_fov, geo_fov[1] / len_geo_fov, geo_fov[2] / len_geo_fov)

def julian(epoch):
    """
    Converts from year/month/day/hour/minute/second/millisecond to Julian days.

    Parameters:
    epoch (tuple): Epoch time as (year, month, day, hour, minute, second, millisecond)

    Returns:
    tuple: Integer part of the Julian day (j) and fractional part of the Julian day (d)
    """
    months = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
    epoch_year = 1988
    epoch_j = 2447160
    epoch_d = 0.5

    year = epoch.year - epoch_year
    leap_year = year // 4
    leap_check = year - 4 * leap_year
    j = 365 * year + months[epoch.month-1] + epoch.day + leap_year + 1
    if leap_check == 0 and epoch.month <= 2:
        j -= 1

    d = (epoch.hour + (epoch.minute + (epoch.second  / 60.0)) / 60.0) / 24.0

    j += epoch_j
    d += epoch_d
    if d > 1:
        d -= 1
        j += 1

    return j, d

def kepler(j, d, elem):
    """
    Computes the altitude, position, and velocity of a satellite at a given time.

    Parameters:
    j (int): Julian day
    d (float): Decimal fraction of a day
    elem (OrbitalElements): Orbital elements

    Returns:
    SatelliteAxes: Object containing satellite's position and velocity vectors
    """
    ERROR_TOL=1e-8
    t = j - elem.jul + d - elem.day
    mean_anom = t * elem.anom_freq + elem.mean_anom

    # Solve Kepler's equation for the eccentric anomaly
    ecc_anom = mean_anom
    old_anom = mean_anom + 1
    while abs(ecc_anom - old_anom) > ERROR_TOL:
        old_anom = ecc_anom
        ecc_anom = mean_anom + elem.ecc * math.sin(old_anom)

    cos_ecc = math.cos(ecc_anom)
    sin_ecc = math.sin(ecc_anom)
    q = math.sqrt((1 - elem.ecc) * (1 + elem.ecc))

    # Compute the coordinates of the satellite
    sat = SatelliteAxes()
    sat.u = (elem.semi_maj_ax * (cos_ecc - elem.ecc),
             elem.semi_maj_ax * sin_ecc * q,
             0.0)
    sat.r = math.sqrt(sum(x ** 2 for x in sat.u))
    sat.u = tuple(x / sat.r for x in sat.u)

    # Calculate the vector V in the plane of T and U, orthogonal to U
    sat.v = (-sat.u[1], sat.u[0], 0.0)

    # Complete the right-handed coordinate system U, V, and W
    sat.w = (0.0, 0.0, 1.0)

    # Calculate the velocity vector
    #print(elem.anom_freq)
    x = elem.anom_freq / (1 - elem.ecc * cos_ecc)
    x *= elem.semi_maj_ax
    #print(x)
    sat.t = (-x * sin_ecc, x * cos_ecc * q, 0.0)
    #print(cos_ecc)
    #print(sat.t)
    sat.s = math.sqrt(sum(x ** 2 for x in sat.t))
    sat.t = tuple(x / sat.s for x in sat.t)

    # Rotate the axes
    perigee = elem.arg_perigee + t * elem.prec_perigee
    cos_rot = math.cos(perigee)
    sin_rot = math.sin(perigee)

    sat.t = (cos_rot * sat.t[0] - sin_rot * sat.t[1],
             sin_rot * sat.t[0] + cos_rot * sat.t[1],
             sat.t[2])
    sat.u = (cos_rot * sat.u[0] - sin_rot * sat.u[1],
             sin_rot * sat.u[0] + cos_rot * sat.u[1],
             sat.u[2])
    sat.v = (cos_rot * sat.v[0] - sin_rot * sat.v[1],
             sin_rot * sat.v[0] + cos_rot * sat.v[1],
             sat.v[2])
    sat.w = (cos_rot * sat.w[0] - sin_rot * sat.w[1],
             sin_rot * sat.w[0] + cos_rot * sat.w[1],
             sat.w[2])

    cos_rot = math.cos(elem.inc)
    sin_rot = math.sin(elem.inc)

    sat.t = (sat.t[0],
             cos_rot * sat.t[1] - sin_rot * sat.t[2],
             sin_rot * sat.t[1] + cos_rot * sat.t[2])
    sat.u = (sat.u[0],
             cos_rot * sat.u[1] - sin_rot * sat.u[2],
             sin_rot * sat.u[1] + cos_rot * sat.u[2])
    sat.v = (sat.v[0],
             cos_rot * sat.v[1] - sin_rot * sat.v[2],
             sin_rot * sat.v[1] + cos_rot * sat.v[2])
    sat.w = (sat.w[0],
             cos_rot * sat.w[1] - sin_rot * sat.w[2],
             sin_rot * sat.w[1] + cos_rot * sat.w[2])

    asc_node = elem.asc_node + t * elem.prec_asc_node
    cos_rot = math.cos(asc_node)
    sin_rot = math.sin(asc_node)

    sat.t = (cos_rot * sat.t[0] - sin_rot * sat.t[1],
             sin_rot * sat.t[0] + cos_rot * sat.t[1],
             sat.t[2])
    sat.u = (cos_rot * sat.u[0] - sin_rot * sat.u[1],
             sin_rot * sat.u[0] + cos_rot * sat.u[1],
             sat.u[2])
    sat.v = (cos_rot * sat.v[0] - sin_rot * sat.v[1],
             sin_rot * sat.v[0] + cos_rot * sat.v[1],
             sat.v[2])
    sat.w = (cos_rot * sat.w[0] - sin_rot * sat.w[1],
             sin_rot * sat.w[0] + cos_rot * sat.w[1],
             sat.w[2])

    return sat

def nav_setup():
    """
    Initializes the sidereal and solar epochs, anomalistic frequencies, and the ((a/b)**2 - 1) quantity for the Earth.
    """
    greenwich_elements.jul, greenwich_elements.day = julian(greenwich_elements.epoch)
    sun_elements.jul, sun_elements.day = julian(sun_elements.epoch)

    greenwich_elements.anom_freq = TWO_PI / greenwich_elements.anom_period
    sun_elements.anom_freq = TWO_PI / sun_elements.anom_period

    global GEO_FACTOR
    GEO_FACTOR = (GEO_ECC / (1 - GEO_ECC)) * (GEO_ECC / (1 + GEO_ECC))

def geocentric_coordinates(j, d, lat, lon):
    """
    Computes the geocentric coordinates of a point with specified terrestrial coordinates.

    Args:
        j (int): Julian day
        d (float): Fraction of Julian day
        lat (float): Geodetic latitude of the point (radians)
        lon (float): Longitude of the point (radians)

    Returns:
        tuple: Unit vector from the center of the earth to the point on the surface
    """
    asc_node = greenwich_elements.anom_freq * ((j - greenwich_elements.jul) + (d - greenwich_elements.day))
    phi = lon + asc_node
    lat_tol = 1.0e-10

    if lat > lat_tol:
        theta = atan((1 + GEO_FACTOR) / tan(lat))
    elif lat < -lat_tol:
        theta = atan((1 + GEO_FACTOR) / tan(lat)) + pi
    else:
        theta = pi / 2

    geo_pos = (sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
    return geo_pos

def terrestrial_coordinates(j, d, geo_pos):
    """
    Computes the terrestrial coordinates of a point with specified geocentric coordinates.

    Args:
        j (int): Julian day
        d (float): Fraction of Julian day
        geo_pos (tuple): Unit vector from the center of the earth to the point on the surface

    Returns:
        tuple: Latitude and longitude of the point (radians)
    """
    asc_node = greenwich_elements.anom_freq * ((j - greenwich_elements.jul) + (d - greenwich_elements.day))
    lat = atan2(geo_pos[2] * (1 + GEO_FACTOR), sqrt((1 - geo_pos[2]) * (1 + geo_pos[2])))
    lon = atan2(geo_pos[1], geo_pos[0]) - asc_node

    lon = lon % (2 * pi)
    if lon < 0:
        lon += 2 * pi

    return lat, lon

def right_ascension_greenwich(j, d):
    """
    Computes the right ascension of Greenwich at the time specified by Julian day and fractional day.

    Args:
        j (int): Julian day
        d (float): Fraction of Julian day

    Returns:
        float: Right ascension of Greenwich
    """
    return greenwich_elements.anom_freq * ((j - greenwich_elements.jul) + (d - greenwich_elements.day))

def solar_position(j, d, new_station, station_lat, station_lon):
    """
    Computes the solar zenith and azimuth from universal time and local latitude and longitude.

    Args:
        j (int): Julian day
        d (float): Decimal fraction of a day
        new_station (bool): Flag indicating if station latitude and longitude have changed
        station_lat (float): Station geodetic latitude (radians)
        station_lon (float): Station longitude (radians)

    Returns:
        tuple: Solar zenith angle (radians), solar azimuth angle (radians), and solar position (SatelliteAxes)
    """
    if new_station:
        new_station = False
        co_lat = pi / 2 - atan(tan(station_lat) * (1 - GEO_ECC) * (1 + GEO_ECC))
        cos_co_lat = cos(co_lat)
        sin_co_lat = sin(co_lat)

    local_sun = kepler(j, d, sun_elements)

    asc_node = right_ascension_greenwich(j, d)
    azimuth = station_lon + asc_node

    cos_azimuth = cos(azimuth)
    sin_azimuth = sin(azimuth)
    pos = (sin_co_lat * cos_azimuth, sin_co_lat * sin_azimuth, cos_co_lat)
    north = (-cos_co_lat * cos_azimuth, -cos_co_lat * sin_azimuth, sin_co_lat)
    east = (-sin_azimuth, cos_azimuth, 0)

    sun_zenith = dot_product(pos, local_sun.u)
    sun_bearing = (local_sun.u[0] - sun_zenith * pos[0],
                   local_sun.u[1] - sun_zenith * pos[1],
                   local_sun.u[2] - sun_zenith * pos[2])
    sun_bearing = [x / sqrt(sum(x ** 2 for x in sun_bearing)) for x in sun_bearing]

    sun_azimuth = dot_product(sun_bearing, north)

    sun_zenith = acos(sun_zenith)
    sun_azimuth = acos(sun_azimuth)
    if dot_product(sun_bearing, east) < 0:
        sun_azimuth = 2 * pi - sun_azimuth

    return sun_zenith, sun_azimuth, local_sun

def satellite_position(j, d, new_station, station_lat, station_lon, sat_elements):
    if new_station:
        new_station = False
        geodetic_lat = station_lat
        geodetic_co_lat = PI / 2 - geodetic_lat
        cos_geodetic_co_lat = math.cos(geodetic_co_lat)
        sin_geodetic_co_lat = math.sin(geodetic_co_lat)

        geocentric_lat = math.atan(math.tan(geodetic_lat) * (1 - GEO_ECC) / (1 + GEO_ECC))
        geocentric_co_lat = PI / 2 - geocentric_lat
        cos_geocentric_co_lat = math.cos(geocentric_co_lat)
        sin_geocentric_co_lat = math.sin(geocentric_co_lat)

    local_sat = kepler(j, d, sat_elements)

    asc_node = greenwich_elements.anom_freq * ((j - greenwich_elements.jul) + (d - greenwich_elements.day))

    fov_azimuth = station_lon + asc_node
    cos_fov_azimuth = math.cos(fov_azimuth)
    sin_fov_azimuth = math.sin(fov_azimuth)

    geo_fov = (sin_geocentric_co_lat * cos_fov_azimuth,
               sin_geocentric_co_lat * sin_fov_azimuth,
               cos_geocentric_co_lat)

    len_geo_fov = geo_radius(geo_fov)

    fov_sat = [local_sat.r * u - len_geo_fov * geo for u, geo in zip(local_sat.u, geo_fov)]

    len_fov_sat = math.sqrt(sum(x ** 2 for x in fov_sat))
    fov_sat = [x / len_fov_sat for x in fov_sat]

    normal = (sin_geodetic_co_lat * cos_fov_azimuth,
              sin_geodetic_co_lat * sin_fov_azimuth,
              cos_geodetic_co_lat)

    north = (-cos_geodetic_co_lat * cos_fov_azimuth,
             -cos_geodetic_co_lat * sin_fov_azimuth,
             sin_geodetic_co_lat)

    east = (-sin_fov_azimuth, cos_fov_azimuth, 0)

    cos_sat_zenith = sum(n * f for n, f in zip(normal, fov_sat))

    sat_bearing = [f - cos_sat_zenith * n for f, n in zip(fov_sat, normal)]
    len_sat_bearing = math.sqrt(sum(x ** 2 for x in sat_bearing))
    sat_bearing = [x / len_sat_bearing for x in sat_bearing]

    cos_sat_azimuth = sum(n * s for n, s in zip(north, sat_bearing))

    sat_zenith = math.acos(cos_sat_zenith)
    sat_azimuth = math.acos(cos_sat_azimuth)
    if sum(e * s for e, s in zip(east, sat_bearing)) < 0:
        sat_azimuth = TWO_PI - sat_azimuth

    return sat_zenith, sat_azimuth, new_station

def element_unit_conversion(elem):
    elem.mean_anom *= RADIANS_PER_DEGREE
    elem.inc *= RADIANS_PER_DEGREE
    elem.arg_perigee *= RADIANS_PER_DEGREE
    elem.asc_node *= RADIANS_PER_DEGREE

    elem.period = (TWO_PI / GEO_GRAV) * elem.semi_maj_ax ** 1.5
    elem.period /= SECONDS_PER_DAY
    elem.freq = TWO_PI / elem.period

    if elem.anom_period == 0:
        elem.anom_period = 0
        elem.anom_freq = 0
    else:
        elem.anom_period /= SECONDS_PER_DAY
        elem.anom_freq = TWO_PI / elem.anom_period

    elem.prec_perigee *= RADIANS_PER_DEGREE
    elem.prec_asc_node *= RADIANS_PER_DEGREE

    # Calculate Julian day and fractional day at epoch
    elem.jul, elem.day = julian(elem.epoch)

    return elem

def perturb_elements(elem, order):
    if order == 0:
        elem.anom_freq = elem.freq
        elem.anom_period = elem.period
        elem.prec_perigee = 0
        elem.prec_asc_node = 0
    elif order == 2:
        sin_inc = math.sin(elem.inc)
        cos_inc = math.cos(elem.inc)
        x = (1 - elem.ecc) * (1 + elem.ecc)
        p = x * elem.semi_maj_ax / GEO_SEMIMAJAX
        q = 1.5 * J_2 / p ** 2
        r = sin_inc ** 2

        elem.anom_freq = elem.freq * (1 + q * math.sqrt(x) * (1 - 1.5 * r))
        elem.anom_period = TWO_PI / elem.anom_freq
        elem.prec_perigee = elem.anom_freq * q * (2 - 2.5 * r)
        elem.prec_asc_node = -elem.anom_freq * q * cos_inc
    else:
        print("Illegal order in Perturb_Elements!")
        print(f"Order = {order}")

    return elem

def cross_product(x, y):
    """
    Evaluates the cross product of two 3D vectors X and Y.

    cross_product = X x Y

    Args:
        x (tuple[float, float, float]): Vector X
        y (tuple[float, float, float]): Vector Y

    Returns:
        tuple[float, float, float]: Cross product of X and Y
    """
    z = (
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    )
    return z


def triple_product(x, y, z):
    """
    Evaluates the triple product of three 3D vectors X, Y, and Z.

    triple_product = X . (Y x Z)

    Args:
        x (tuple[float, float, float]): Vector X
        y (tuple[float, float, float]): Vector Y
        z (tuple[float, float, float]): Vector Z

    Returns:
        float: Triple product of X, Y, and Z
    """
    triple_product = (
        x[0] * (y[1] * z[2] - y[2] * z[1])
        + x[1] * (y[2] * z[0] - y[0] * z[2])
        + x[2] * (y[0] * z[1] - y[1] * z[0])
    )
    return triple_product
