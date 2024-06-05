import math

# Constants
DEGREES_PER_RADIAN = 180.0 / math.pi
RADIANS_PER_DEGREE = math.pi / 180.0
TWO_PI = 2.0 * math.pi

class SolarPositionData:
    """Solar position computed by Solar_Position_Almanac"""
    def __init__(self):
        self.GeoSun = [0.0, 0.0, 0.0]  # Unit vector from Earth to the sun
        self.DSunEarth = 0.0  # Sun-Earth distance in Astronomical Units
        self.RightAscension = 0.0  # Solar right ascension
        self.Declination = 0.0  # Solar declination
        self.EqnTime = 0.0  # Apparent minus mean time (minutes)

def write_xminsec(string, angle, ident=None):
    """Writes out an angle (time) as degrees (hours), minutes, seconds"""
    if ident:
        angle_deg = angle
    else:
        angle_deg = angle * DEGREES_PER_RADIAN
    degrees = int(angle_deg)
    angle_min = (angle_deg - degrees) * 60
    minutes = int(angle_min)
    seconds = (angle_min - minutes) * 60
    print(f"{string:30}{degrees:4}{minutes:3}{seconds:6.2f}")

def apparent_zenith(true_zenith):
    """Calculates the apparent solar zenith distance from the true solar zenith distance (in degrees)"""
    c = math.cos(true_zenith * RADIANS_PER_DEGREE)
    if c > 0:
        d = -0.0471210419 + 1.0 / (0.955000 + 20.267 * c)
        c = c + 0.008300 * d
        apparent_zenith = math.acos(c) * DEGREES_PER_RADIAN
    else:
        apparent_zenith = 90.0
    return apparent_zenith

def solar_position_almanac(jul_day, tim_day):
    """Computes the solar position vector in equatorial coordinates"""
    data = SolarPositionData()
    day_2000 = jul_day + tim_day - 2451545

    # Calculate mean longitude of sun
    mean_lon = (280.472 + 0.9856474 * day_2000) * RADIANS_PER_DEGREE
    mean_lon = put_in_range(mean_lon, TWO_PI)

    # Calculate mean anomaly
    mean_anom = (357.528 + 0.9856003 * day_2000) * RADIANS_PER_DEGREE
    mean_anom = put_in_range(mean_anom, TWO_PI)
    g = mean_anom

    # Calculate ecliptic longitude of sun
    ecliptic_lon = mean_lon + (1.915 * math.sin(g) + 0.020 * math.sin(2 * g)) * RADIANS_PER_DEGREE
    ecliptic_lon = put_in_range(ecliptic_lon, TWO_PI)

    # Calculate obliquity of the ecliptic
    obliquity_ecliptic = (23.439 - 0.0000004 * day_2000) * RADIANS_PER_DEGREE

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
    if x < -math.pi:
        x = x + TWO_PI
    elif x > math.pi:
        x = x - TWO_PI
    data.EqnTime = x * DEGREES_PER_RADIAN * 4.0

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

def atan_quad(x, angle):
    """Computes the inverse tangent of x corrected to lie in the same quadrant as angle"""
    i = quadrant(angle)
    atan_quad = math.atan(x)
    if i == 2 or i == 3:
        atan_quad = atan_quad + math.pi
    elif i == 4:
        atan_quad = atan_quad + TWO_PI
    if abs(atan_quad - angle) >= math.pi / 4:
        print("Failure in Atan_Quad!")
    return atan_quad

def quadrant(x):
    """Returns the quadrant containing the input angle (in radians)"""
    if x < 0.0 or x > TWO_PI:
        print("Trouble in Quadrant, argument out of range!")
    if 0.0 <= x < 0.5 * math.pi:
        return 1
    elif 0.5 * math.pi <= x < math.pi:
        return 2
    elif math.pi <= x < 1.5 * math.pi:
        return 3
    elif 1.5 * math.pi <= x < 2.0 * math.pi:
        return 4

def normalize_vector(vector):
    """Normalizes a vector to unit length"""
    length = math.sqrt(sum(x ** 2 for x in vector))
    if length > 0:
        return [x / length for x in vector]
    else:
        return vector