import datetime as dt
import math

pi = math.pi
rad2deg = 180.0 / pi
deg2rad = pi /180.0
earthrad = 6370.0 #km
gravity = 9.80665  # ms^{-2}, standard value

# Physical constants for the earth
#GEO_MASS = 5.9733e+24  # kg
GEO_SEMI_MAJ_AX = 6378.388e+03  # m
GEO_SEMI_MIN_AX = 6356.912e+03  # m
GEO_ECC = 8.199179e-02
#GEO_FLATTENING = 3.367e-03
#GEO_GRAV = 1.9964985e+07  # m^{1.5}/s  

def geo_radius(geo_pos):
    return GEO_SEMI_MAJ_AX / math.sqrt(1 + GEO_FACTOR * geo_pos[2] ** 2)

class SolarPositionData:
    """Solar position computed by Solar_Position_Almanac"""
    def __init__(self):
        self.GeoSun = [0.0, 0.0, 0.0]  # Unit vector from Earth to the sun
        self.DSunEarth = 0.0  # Sun-Earth distance in Astronomical Units
        self.RightAscension = 0.0  # Solar right ascension
        self.Declination = 0.0  # Solar declination
        self.EqnTime = 0.0  # Apparent minus mean time (minutes)

def put_in_range(x, y):
    """Ensures that x ends up in the range (0, y) by incrementing in units of y"""
    if x < 0.0:
        while x < 0.0:
            x = x + y
    if x >= y:
        while x >= y:
            x = x - y
    return x

def atan_quad(x, angle):
    """Computes the inverse tangent of x corrected to lie in the same quadrant as angle"""
    i = quadrant(angle)
    atan_quad = math.atan(x)
    if i == 2 or i == 3:
        atan_quad = atan_quad + pi
    elif i == 4:
        atan_quad = atan_quad + 2*pi
    if abs(atan_quad - angle) >= (pi / 4):
        print("Failure in Atan_Quad!")
    return atan_quad

def quadrant(x):
    """Returns the quadrant containing the input angle (in radians)"""
    if x < 0.0 or x > (2*pi):
        print("Trouble in Quadrant, argument out of range!")
    if 0.0 <= x < 0.5 * pi:
        return 1
    elif 0.5 * pi <= x < pi:
        return 2
    elif pi <= x < 1.5 * pi:
        return 3
    elif 1.5 * pi <= x < (2.0 * pi):
        return 4

def normalize_vector(vector):
    """Normalizes a vector to unit length"""
    length = math.sqrt(sum(x ** 2 for x in vector))
    if length > 0:
        return [x / length for x in vector]
    else:
        return vector
def solar_position_almanac(jul_day, tim_day):
    """Computes the solar position vector in equatorial coordinates"""
    data = SolarPositionData()
    day_2000 = jul_day + tim_day - 2451545
    #print('day2000 = {}'.format(day_2000))
    # Calculate mean longitude of sun
    mean_lon = (280.472 + 0.9856474 * day_2000) * deg2rad
    #print('mean lon1 ',mean_lon)
    mean_lon = put_in_range(mean_lon, 2*pi)
    #print('mean lon2 ',mean_lon)
    # Calculate mean anomaly
    mean_anom = (357.528 + 0.9856003 * day_2000) * deg2rad
    #print('mean_anom = ',mean_anom)
    mean_anom = put_in_range(mean_anom, 2*pi)
    g = mean_anom
    #print('g = {}'.format(g))
    # Calculate ecliptic longitude of sun
    ecliptic_lon = mean_lon + (1.915 * math.sin(g) + 0.020 * math.sin(2 * g)) * deg2rad
    #print('ecl lon1',ecliptic_lon)
    ecliptic_lon = put_in_range(ecliptic_lon, 2*pi)
    #print('ecl lon2',ecliptic_lon)

    # Calculate obliquity of the ecliptic
    obliquity_ecliptic = (23.439 - 0.0000004 * day_2000) * deg2rad

    # Calculate right ascension
    x = math.cos(obliquity_ecliptic) * math.tan(ecliptic_lon)
    data.RightAscension = atan_quad(x, ecliptic_lon)

    # Declination
    x = math.sin(obliquity_ecliptic) * math.sin(ecliptic_lon)
    data.Declination = math.asin(x)

    # Sun Earth distance
    data.DSunEarth = 1.00014 - 0.01671 * math.cos(g) - 0.00014 * math.cos(2 * g)

    # Solar position vector
    data.GeoSun[0] = data.DSunEarth * math.cos(ecliptic_lon)
    data.GeoSun[1] = data.DSunEarth * math.cos(obliquity_ecliptic) * math.sin(ecliptic_lon)
    data.GeoSun[2] = data.DSunEarth * math.sin(obliquity_ecliptic) * math.sin(ecliptic_lon)

    # Normalize the position vector
    r = normalize_vector(data.GeoSun)

    # Equation of time (apparent solar time - mean solar time) (minutes)
    x = mean_lon - data.RightAscension
    if x < -pi:
        x = x + 2*pi
    elif x > pi:
        x = x - 2*pi
    data.EqnTime = x * rad2deg * 4.0

    return data

def put_in_range(x, y):
    """Ensures that x ends up in the range (0, y) by incrementing in units of y"""
    if x < 0.0:
        while x < 0.0:
            x = x + y
    if x >= y:
        while x >= y:
            x = x - y
    return x

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

# Define default orbital elements
default_elements = OrbitalElements()
greenwich_elements = OrbitalElements(
    semi_maj_ax=1.0, period=0.997269566340, anom_period=0.997269566340,
    #epoch=(1988, 1, 0, 17, 21, 35, 354)
    epoch=dt.datetime(1988,1,1,17,21,35,354)
)

sun_elements = OrbitalElements(
    id_code=1, semi_maj_ax=1.496e11, ecc=0.167133900e-01, mean_anom=0.621622191e+01,
    inc=0.409120047e+00, arg_perigee=0.493460338e+01, period=0.365259635e+03,
    anom_period=0.365259635e+03, prec_perigee=0.000000822e+00, 
    #epoch=(1988, 1, 0, 0, 0, 0, 0)
    epoch=dt.datetime(1988,1,1,0,0,0,0) 
)
greenwich_elements.jul, greenwich_elements.day = julian(greenwich_elements.epoch)
sun_elements.jul, sun_elements.day = julian(sun_elements.epoch)

greenwich_elements.anom_freq = 2*pi / greenwich_elements.anom_period
sun_elements.anom_freq = 2*pi / sun_elements.anom_period

GEO_FACTOR = (GEO_ECC / (1 - GEO_ECC)) * (GEO_ECC / (1 + GEO_ECC))




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


def satellite_position(j, d, station_lat, station_lon, sat_elements):

    geodetic_lat = station_lat
    geodetic_co_lat = pi / 2 - geodetic_lat
    cos_geodetic_co_lat = math.cos(geodetic_co_lat)
    sin_geodetic_co_lat = math.sin(geodetic_co_lat)

    geocentric_lat = math.atan(math.tan(geodetic_lat) * (1 - GEO_ECC) / (1 + GEO_ECC))
    geocentric_co_lat = pi / 2 - geocentric_lat
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
        sat_azimuth = 2*pi - sat_azimuth

    return sat_zenith, sat_azimuth


def apparent_zenith(true_zenith):
    """Calculates the apparent solar zenith distance from the true solar zenith distance (in degrees)"""
    c = math.cos(true_zenith*deg2rad)
    if c > 0:
        d = -0.0471210419 + 1.0 / (0.955000 + 20.267 * c)
        c = c + 0.008300 * d
        apparent_zenith = math.acos(c) * rad2deg
    else:
        apparent_zenith = 90.0
    return apparent_zenith