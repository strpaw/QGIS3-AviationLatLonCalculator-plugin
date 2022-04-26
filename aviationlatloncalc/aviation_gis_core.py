import math
import re

""" Constants """

# Parameters of WGS84 ellipsoid
WGS84_A = 6378137.0  # semi-major axis of the WGS84 ellipsoid in m
WGS84_B = 6356752.314245  # semi-minor axis of the WGS84 ellipsoid in m
WGS84_F = 1 / 298.257223563  # flattening of the WGS84 ellipsoid

# Units of measure
UOM_M = 'M'  # meters
UOM_KM = 'KM'  # kilometers
UOM_NM = 'NM'  # nautical miles
UOM_FEET = 'FEET'  # feet
UOM_SM = 'SM'  # statue miles

UOM_LIST = [UOM_M, UOM_KM, UOM_NM, UOM_FEET, UOM_SM]

# Geometry types

GEOM_POINT = 'POINT'
GEOM_LINE = 'LINE'
GEOM_POLYGON = 'POLYGON'

C_LAT = 'C_LAT'
C_LON = 'C_LON'

ANGLE_PREFIX = ['-', '+', 'N', 'S', 'E', 'W']
ANGLE_SUFFIX = ['N', 'S', 'E', 'W']
ANGLE_POSITIVE_SIGN = ['+', 'N', 'E']
ANGLE_NEGATIVE_SIGN = ['-', 'S', 'W']

ANGLE_POSITIVE = 'COORD_POSITIVE'
ANGLE_NEGATIVE = 'COORD_NEGATIVE'


S_SPACE = ' '
S_HYPHEN = '-'
S_DEG_WORD = 'DEG'
S_DEG_LETTER = 'D'
S_MIN_WORD = 'MIN'
S_MIN_LETTER = 'M'
S_SEC_WORD = 'SEC'
S_ALL = [S_SPACE, S_HYPHEN, S_DEG_WORD, S_DEG_LETTER, S_MIN_WORD, S_MIN_LETTER, S_SEC_WORD]


# DMS, DM compacted constant formats
F_DMS_COMP = 'F_DMS_COMP'  # DMS compacted
F_DM_COMP = 'F_DM_COMP'  # DM compacted

""" Ellipsoid calculations """

def vincenty_direct_solution(begin_lat, begin_lon, begin_azimuth, distance, a, b, f):
    """ Computes the latitude and longitude of the second point based on latitude, longitude,
    of the first point and distance and azimuth from first point to second point.
    Uses the algorithm by Thaddeus Vincenty for direct geodetic problem.
    For more information refer to: http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
    :param begin_lon: float, longitude of the first point; decimal degrees
    :param begin_lat: float, latitude of the first point; decimal degrees
    :param begin_azimuth: float, azimuth from first point to second point; decimal degrees
    :param distance: float, distance from first point to second point; meters
    :param a: float, semi-major axis of ellipsoid; meters
    :param b: float, semi-minor axis of ellipsoid; meters
    :param f: float, flattening of ellipsoid
    :return lat2_dd, lon2_dd: float, float latitude and longitude of the second point, decimal degrees
    """
    # Convert latitude, longitude, azimuth of the initial point to radians
    lat1 = math.radians(begin_lat)
    lon1 = math.radians(begin_lon)
    alfa1 = math.radians(begin_azimuth)
    sin_alfa1 = math.sin(alfa1)
    cos_alfa1 = math.cos(alfa1)

    # U1 - reduced latitude
    tan_u1 = (1 - f) * math.tan(lat1)
    cos_u1 = 1 / math.sqrt(1 + tan_u1 * tan_u1)
    sin_u1 = tan_u1 * cos_u1

    # sigma1 - angular distance on the sphere from the equator to initial point
    sigma1 = math.atan2(tan_u1, math.cos(alfa1))

    # sin_alfa - azimuth of the geodesic at the equator
    sin_alfa = cos_u1 * sin_alfa1
    cos_sq_alfa = 1 - sin_alfa * sin_alfa
    u_sq = cos_sq_alfa * (a * a - b * b) / (b * b)
    A = 1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
    B = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))

    sigma = distance / (b * A)
    sigmap = 1
    sin_sigma, cos_sigma, cos2sigma_m = None, None, None

    while math.fabs(sigma - sigmap) > 1e-12:
        cos2sigma_m = math.cos(2 * sigma1 + sigma)
        sin_sigma = math.sin(sigma)
        cos_sigma = math.cos(sigma)
        d_sigma = B * sin_sigma * (cos2sigma_m + B / 4 * (
                    cos_sigma * (-1 + 2 * cos2sigma_m * cos2sigma_m) - B / 6 * cos2sigma_m * (
                        -3 + 4 * sin_sigma * sin_sigma) * (-3 + 4 * cos2sigma_m * cos2sigma_m)))
        sigmap = sigma
        sigma = distance / (b * A) + d_sigma

    var_aux = sin_u1 * sin_sigma - cos_u1 * cos_sigma * cos_alfa1  # Auxiliary variable

    # Latitude of the end point in radians
    lat2 = math.atan2(sin_u1 * cos_sigma + cos_u1 * sin_sigma * cos_alfa1,
                      (1 - f) * math.sqrt(sin_alfa * sin_alfa + var_aux * var_aux))

    lamb = math.atan2(sin_sigma * sin_alfa1, cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos_alfa1)
    C = f / 16 * cos_sq_alfa * (4 + f * (4 - 3 * cos_sq_alfa))
    L = lamb - (1 - C) * f * sin_alfa * (
                sigma + C * sin_sigma * (cos2sigma_m + C * cos_sigma * (-1 + 2 * cos2sigma_m * cos2sigma_m)))
    # Longitude of the second point in radians
    lon2 = (lon1 + L + 3 * math.pi) % (2 * math.pi) - math.pi

    # Convert to decimal degrees
    lat2_dd = math.degrees(lat2)
    lon2_dd = math.degrees(lon2)

    return lat2_dd, lon2_dd


def vincenty_reverse_solution(lat1_dd, lon1_dd, lat2_dd, lon2_dd, a, b, f):
    """ Computes distance and bearing between two points
    :param lat1_dd: float, latitude in decimal degrees
    :param lon1_dd: float, longitude in decimal degrees
    :param lat2_dd: float, latitude in decimal degrees
    :param lon2_dd: float, longitude in decimal degrees
    :param a: float, semi-major axis of ellipsoid; meters
    :param b: float, semi-minor axis of ellipsoid; meters
    :param f: float, flattening of ellipsoid
    :return: alfa1: float: bearing from point 1 to point 2
             distance: distance between point 1 and 2
    """
    lat1 = math.radians(lat1_dd)
    lon1 = math.radians(lon1_dd)
    lat2 = math.radians(lat2_dd)
    lon2 = math.radians(lon2_dd)

    L = lon2 - lon1

    tan_u1 = (1 - f) * math.tan(lat1)
    cos_u1 = 1 / math.sqrt((1 + tan_u1 * tan_u1))
    sin_u1 = tan_u1 * cos_u1
    tan_u2 = (1 - f) * math.tan(lat2)
    cos_u2 = 1 / math.sqrt((1 + tan_u2 * tan_u2))
    sin_u2 = tan_u2 * cos_u2

    lamb = L
    lamb_p = L
    iterations = 0

    cos_sq_alfa, sin_sigma, cos2sigma__m, sin_lamb, cos_lamb, cos_sigma = None, None, None, None, None, None
    sigma = None

    while math.fabs(lamb - lamb_p) > 1e-12 or iterations < 100:
        iterations += 1
        sin_lamb = math.sin(lamb)
        cos_lamb = math.cos(lamb)
        sin_sq_sigma = (cos_u2 * sin_lamb) * (cos_u2 * sin_lamb) + (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lamb) *\
                       (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lamb)
        sin_sigma = math.sqrt(sin_sq_sigma)
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lamb
        sigma = math.atan2(sin_sigma, cos_sigma)
        sin_alfa = cos_u1 * cos_u2 * sin_lamb / sin_sigma
        cos_sq_alfa = 1 - sin_alfa * sin_alfa
        cos2sigma_m = cos_sigma - 2 * sin_u1 * sin_u2 / cos_sq_alfa
        C = f / 16 * cos_sq_alfa * (4 + f * (4 - 3 * cos_sq_alfa))
        lamb_p = lamb
        lamb = L + (1 - C) * f * sin_alfa * (
                    sigma + C * sin_sigma * (cos2sigma_m + C * cos_sigma * (-1 + 2 * cos2sigma_m * cos2sigma_m)))
        iterations += 1

    u_sq = cos_sq_alfa * (a * a - b * b) / (b * b)
    A = 1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
    B = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))
    delta_sigma = B * sin_sigma * (cos2sigma_m + B / 4 * (
                cos_sigma * (-1 + 2 * cos2sigma_m * cos2sigma_m) - B / 6 * cos2sigma_m * (
                    -3 + 4 * sin_sigma * sin_sigma) * (-3 + 4 * cos2sigma_m * cos2sigma_m)))

    distance = b * A * (sigma - delta_sigma)
    alfa1 = math.atan2(cos_u2 * sin_lamb, cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lamb)
    alfa1 = (alfa1 + 2 * math.pi) % (2 * math.pi)
    alfa1 = math.degrees(alfa1)

    return alfa1, distance


def dist_azm_orth_offset2latlon(ref_lat, ref_lon, ref_azm, distance_m, offset_m, offset_side):
    """ Calculates latitude and longitude of the second point base don latitude, longitude of the first point, azimuth,
    distance and orthogonal offset
    Example: distance 1500 m, azimuth 45 degrees and offset 500 meter left
    :param ref_lat: float, reference point latitude
    :param ref_lon: float, reference point longitude
    :param ref_azm: float, azimuth from reference point to intermediate point
    :param distance_m: float, distance in meters
    :param offset_m: float, offset in meters
    :param offset_side: indicate offset side, 'LEFT' for left, 'RIGHT' for right
    :return lat2_dd, lon2_dd: float, second point latitude, longitude
    """
    # Calculate azimuth from intermediate point to second point
    if offset_side == 'LEFT':
        offset_azm = ref_azm - 90
    elif offset_side == 'RIGHT':
        offset_azm = ref_azm + 90
    # Normalize azm to [0,360] degrees
    if offset_azm < 0:
        offset_azm += 360
    elif offset_azm > 360:
        offset_azm -= 360

    # Calculate intermediate point latitude, longitude
    inter_lat_dd, inter_lon_dd = vincenty_direct_solution(ref_lat, ref_lon, ref_azm, distance_m, WGS84_A, WGS84_B,
                                                          WGS84_F)

    # Calculate second point latitude, longitude, as reference point use intermediate point
    lat2_dd, lon2_dd = vincenty_direct_solution(inter_lat_dd, inter_lon_dd, offset_azm, offset_m, WGS84_A, WGS84_B,
                                                WGS84_F)

    return lat2_dd, lon2_dd


def dist_azm_orth_offset2latlon2(ref_lat, ref_lon, x_azm, y_azm, x_m, y_m):
    """ Calculates latitude and longitude of the second point based on latitude, longitude of the reference point,
    azimuth, distance and orthogonal offset.
    :param ref_lat: float, reference point latitude
    :param ref_lon: float, reference point longitude
    :param x_azm: float, X axis azimuth
    :param y_azm: float, Y axis azimuth
    :param x_m: float, x coordinate in meters
    :param y_m: float, y coordinate in meters
    :return lat2_dd, lon2_dd: float, second point latitude, longitude
    """
    # Calculate intermediate point latitude, longitude
    inter_lat_dd, inter_lon_dd = vincenty_direct_solution(ref_lat, ref_lon, x_azm, x_m, WGS84_A, WGS84_B,
                                                          WGS84_F)

    # Calculate second point latitude, longitude, as reference point use intermediate point
    lat2_dd, lon2_dd = vincenty_direct_solution(inter_lat_dd, inter_lon_dd, y_azm, y_m, WGS84_A, WGS84_B,
                                                WGS84_F)

    return lat2_dd, lon2_dd


""" Base """


class AviationBaseClass:
    """ Aviation base class for storing and manipulating data with aviation content """

    def __init__(self):
        self._is_valid = None
        self._err_msg = ''

    @staticmethod
    def normalize_src_input(src_input):
        """ Normalizes source (input)  value for further processing
        :param src_input: str, input angle string to normalize
        :return: norm_angle: str, normalized angle string
        """
        norm_input = str(src_input).strip()
        norm_input = norm_input.replace(',', '.')
        norm_input = norm_input.upper()
        return norm_input

    @property
    def is_valid(self):
        return self._is_valid

    @is_valid.setter
    def is_valid(self, value):
        self._is_valid = value

    @property
    def err_msg(self):
        return self._err_msg

    @err_msg.setter
    def err_msg(self, value):
        self._err_msg = value


""" Distance """


class Distance(AviationBaseClass):
    def __init__(self, dist_src, dist_uom, checked_value='Distance', allow_negative=False):
        AviationBaseClass.__init__(self)
        self.dist_src = dist_src  # Source distance
        self.dist_uom = dist_uom  # Unit of measure in distance is expressed which
        self.allow_negative = allow_negative  # Indicates if negative values are allowed (e. g. Cartesian)
        self._dist_float = None  # Distance float value
        self.checked_value = checked_value
        self.check_distance()

    @staticmethod
    def to_meters(d, from_unit):
        """ Converts distance given in specified unit to distance in meters
        :param d: float, distance in unit specified by parameter from_unit
        :param from_unit: constant unit of measure, unit of measure parameter d_unit
        :return float, distance in unit specified by parameter to_unit
        """
        if from_unit == UOM_M:
            return d
        elif from_unit == UOM_KM:
            return d * 1000
        elif from_unit == UOM_NM:
            return d * 1852
        elif from_unit == UOM_FEET:
            return d * 0.3048
        elif from_unit == UOM_SM:
            return d * 1609.344
        else:
            return None

    @staticmethod
    def from_meters(d, to_unit):
        """ Converts distance given in meters to distance in specified unit
        :param d: float, distance in meters
        :param to_unit: constant unit of measurement
        :return float, distance in unit specified by parameter to_unit
        """
        if to_unit == UOM_M:
            return d
        elif to_unit == UOM_KM:
            return d / 1000
        elif to_unit == UOM_NM:
            return d / 1852
        elif to_unit == UOM_FEET:
            return d / 0.3048
        elif to_unit == UOM_SM:
            return d / 1609.344
        else:
            return None

    def convert_distance(self, d, from_unit, to_unit):
        """ Convert distance between various units
        :param d: float, distance in units specified by parameter from_unit
        :param from_unit: constant measure of units
        :param to_unit: constant measure of unit
        :return float, distance in units specified by parameter to_unit
        """
        if from_unit == to_unit:
            return d
        else:
            d_m = self.to_meters(d, from_unit)  # Convert to meters
            return self.from_meters(d_m, to_unit)  # Convert from meters

    def get_meters(self):
        """ Returns source distance in meters """
        if self.is_valid is True:
            return self.to_meters(self.dist_float, self.dist_uom)

    def to_unit(self, to_unit):
        """ Returns source distance in UOM given in to_unit """
        if self.is_valid is True:
            return self.convert_distance(self.dist_float, self.dist_uom, to_unit)

    def check_distance(self):
        """ Distance validation. Uses float() function to check if distance value is a number """
        if self.dist_src == '':
            self.is_valid = False
            self.err_msg = 'Enter ' + self.checked_value + '\n'
        else:
            try:
                dist_norm = float(self.normalize_src_input(self.dist_src))
                if self.allow_negative is False:
                    if dist_norm < 0:  # Check if distance is less than 0
                        self.is_valid = False
                        self.err_msg = self.checked_value + ' error.\n'

                    else:
                        self.is_valid = True
                        self.dist_float = dist_norm
                else:
                    self.is_valid = True
                    self.dist_float = dist_norm
            except ValueError:
                self.is_valid = False
                self.err_msg = self.checked_value + ' error.\n'

    def get_distance_str_info_data(self):
        """ Returns string with information: distance string value and UOM
        Useful when we want to add source information in output """
        dist_str = '{} {}'.format(self.dist_src, self.dist_uom)
        return dist_str

    @property
    def dist_float(self):
        return self._dist_float

    @dist_float.setter
    def dist_float(self, value):
        self._dist_float = value



""" Angle """

# Regular expression patterns for DMS and DM formats
coord_lat_comp_regex = {F_DMS_COMP: re.compile(r'''(?P<deg>^\d{2})  # Degrees
                                                   (?P<min>\d{2})  # Minutes
                                                   (?P<sec>\d{2}(\.\d+)?$)  # Seconds 
                                                ''', re.VERBOSE),
                        F_DM_COMP: re.compile(r'''(?P<deg>^\d{2})  # Degrees
                                                  (?P<min>\d{2}(\.\d+)?$)   # Minutes    
                                              ''', re.VERBOSE)}

coord_lon_comp_regex = {F_DMS_COMP: re.compile(r'''(?P<deg>^\d{3})  # Degrees
                                                   (?P<min>\d{2})  # Minutes
                                                   (?P<sec>\d{2}\.\d+$|\d{2}$)  # Seconds 
                                                ''', re.VERBOSE),
                        F_DM_COMP: re.compile(r'''(?P<deg>^\d{3})  # Degrees
                                                  (?P<min>\d{2}\.\d+$|\d{2}$)   # Minutes    
                                              ''', re.VERBOSE)}




class AngleBase(AviationBaseClass):
    def __init__(self):
        AviationBaseClass.__init__(self)

    @staticmethod
    def get_angle_parts(angle_norm):
        """
        :param angle_norm: str, normalized angle
        :return: tuple:
        """
        hem_char = angle_norm[0]
        if hem_char in ANGLE_PREFIX:
            if hem_char in ANGLE_POSITIVE_SIGN:
                return ANGLE_POSITIVE, angle_norm[1:].strip(), hem_char
            elif hem_char in ANGLE_NEGATIVE_SIGN:
                return ANGLE_NEGATIVE, angle_norm[1:].strip(), hem_char
        else:
            hem_char = angle_norm[-1]
            if hem_char in ANGLE_SUFFIX:
                if hem_char in ANGLE_POSITIVE_SIGN:
                    return ANGLE_POSITIVE, angle_norm[:-1].strip(), hem_char
                elif hem_char in ANGLE_NEGATIVE_SIGN:
                    return ANGLE_NEGATIVE, angle_norm[:-1].strip(), hem_char
            else:
                # Begins with digit -> positive sign
                return ANGLE_POSITIVE, angle_norm.strip(), hem_char

    # @staticmethod
    # def check_angle_range(angle_dd, min_value, max_value):
    #     """ Checks if angle is within closed interval <min_value, max_value>
    #     :param angle_dd: float, angle value to check
    #     :param min_value: float, minimum value
    #     :param max_value: float, maximum value
    #     :return: tuple (bool, float) if angle is within the range
    #              tuple (bool, None) if angle is out of range
    #     """
    #     if min_value <= angle_dd <= max_value:
    #         return True, angle_dd
    #     else:
    #         return False, None

    @staticmethod
    def check_angle_dd(angle_norm):
        """ Checks if angle is in DD format.
        :param angle_norm: float: angle to check
        :return: float, vale of angle if angle is integer of float, const NOT_VALID otherwise
        """
        try:
            dd = float(angle_norm)
            return True, dd
        except ValueError:
            return False, None

    @staticmethod
    def parse_compacted_formats(regex_patterns, coord_part):
        """ Converts latitude or longitude in DMSH format into DD format.
        :param regex_patterns: dictionary of regex object, patterns of DMS formats
        :param coord_part: str, angle to check
        :return: flag, bool,
        :return: dd:, float if DMS is valid format, None otherwise
        :return: coord_format: DMS format constant in which input is if input is valid, None otherwise
        """
        result = False
        dd = None
        for pattern in regex_patterns:  # Check if input matches any pattern
            if regex_patterns.get(pattern).match(coord_part):
                if pattern == F_DMS_COMP:
                    groups = regex_patterns.get(pattern).search(coord_part)
                    d = float(groups.group('deg'))
                    m = float(groups.group('min'))
                    s = float(groups.group('sec'))

                    if m < 60 and s < 60:
                        result = True
                        dd = d + m / 60 + s / 3600

                if pattern == F_DM_COMP:
                    groups = regex_patterns.get(pattern).search(coord_part)
                    d = float(groups.group('deg'))
                    m = float(groups.group('min'))

                    if m < 60:
                        result = True
                        dd = d + m / 60

        return result, dd

    @staticmethod
    def parse_separated_formats(coord):
        # Replace separators (delimiters) with blank (space)
        for sep in S_ALL:
            coord = re.sub(sep, ' ', coord)
        # Replace multiple spaces into single spaces
        coord_mod = re.sub('\s+', ' ', coord)
        c_parts = coord_mod.split(' ')
        if len(c_parts) == 3:  # Assume format DMS separated


            try:
                d = int(c_parts[0])
                if d < 0:
                    return False, None
            except ValueError:
                return False, None

            try:
                m = int(c_parts[1])
                if m < 0 or m >= 60:
                    return False, None
            except ValueError:
                return False, None

            try:
                s = float(c_parts[2])
                if s < 0 or s >= 60:
                    return False, None
            except ValueError:
                return False, None

            try:
                dd = float(d) + float(m) / 60 + s / 3600
                return True, dd
            except ValueError:
                return False, None
        elif len(c_parts) == 2:  # Assume format DM separated

            try:
                d = int(c_parts[0])
                if d < 0:
                    return False, None
            except ValueError:
                return False, None

            try:
                m = int(c_parts[1])
                if m < 0 or m >= 60:
                    return False, None
            except ValueError:
                return False, None

            try:
                dd = float(d) + float(m) / 60
                return True, dd
            except ValueError:
                return False, None

        elif len(c_parts) == 1:  # Assume format DD separated
            try:
                d = float(c_parts[0])
                if d < 0:
                    return False, None
                else:
                    return True, d
            except ValueError:
                return False, None

        else:
            return False, None

    @staticmethod
    def dd2dms(angle_dd):
        """ Convert decimal degrees format into decimal, minutes, seconds format"""
        dd = math.fabs(angle_dd)
        deg = int(math.trunc(dd))
        min_dd = round((dd - deg) * 60, 8)
        minutes = int(math.trunc(min_dd))
        sec = round((min_dd - minutes) * 60, 3)

        dms = '{deg:0>2} {minutes:0>2} {sec:05.3f}'.format(deg=deg, minutes=minutes, sec=sec)

        return dms


class CoordinatesPair(AngleBase):
    def __init__(self, lat_src, lon_src):
        AngleBase.__init__(self)
        self.lat_src = lat_src
        self.lon_src = lon_src
        self._lat_dd = None
        self._lon_dd = None
        self.parse_coordinates2dd()

    def parse_coordinates2dd(self):
        lat_valid, lon_valid = False, False

        if self.lat_src == '':  # Blank input
            self.err_msg += 'Enter latitude value!\n'
        else:
            lat_src_norm = self.normalize_src_input(self.lat_src)
            # Get parts of the latitude
            lat_sign, lat_deg_part, lat_hem = self.get_angle_parts(lat_src_norm)
            if lat_hem in ['W', 'E']:
                lat_valid = False
            else:
                # Check DMS, DM compacted formats
                lat_valid, lat_dd = self.parse_compacted_formats(coord_lat_comp_regex, lat_deg_part)
                if lat_valid is False:
                    # Check separated formats
                    lat_valid, lat_dd = self.parse_separated_formats(lat_deg_part)

                if lat_valid is True:
                    if lat_sign is ANGLE_NEGATIVE:
                        lat_dd = -1 * lat_dd
                    lat_valid, self.lat_dd = self.check_angle_range(lat_dd, -90, 90)

            if lat_valid is False:
                self.err_msg += 'Latitude error!\n'

        if self.lon_src == '':  # Blank input
            self.err_msg += 'Enter longitude value!\n'
        else:
            lon_src_norm = self.normalize_src_input(self.lon_src)
            # Get parts of the latitude
            lon_sign, lon_deg_part, lon_hem = self.get_angle_parts(lon_src_norm)
            if lon_hem in ['N', 'S']:
                lon_valid = False
            else:
                # Check DMS, DM compacted formats
                lon_valid, lon_dd = self.parse_compacted_formats(coord_lon_comp_regex, lon_deg_part)
                if lon_valid is False:
                    # Check separated formats
                    lon_valid, lon_dd = self.parse_separated_formats(lon_deg_part)

                if lon_valid is True:
                    if lon_sign is ANGLE_NEGATIVE:
                        lon_dd = -1 * lon_dd
                    lon_valid, self.lon_dd = self.check_angle_range(lon_dd, -180, 180)

            if lon_valid is False:
                self.err_msg += 'Longitude error!\n'

        if lat_valid is False or lon_valid is False:
            self.is_valid = False
        else:
            self.is_valid = True

    @property
    def lat_dd(self):
        return self._lat_dd

    @lat_dd.setter
    def lat_dd(self, value):
        self._lat_dd = value

    @property
    def lon_dd(self):
        return self._lon_dd

    @lon_dd.setter
    def lon_dd(self, value):
        self._lon_dd = value


class MagVar(AngleBase):
    def __init__(self, mag_var_src):
        AngleBase.__init__(self)
        self.mag_var_src = mag_var_src
        self._mag_var_dd = None

        self.parse_mag_var2dd()

    def parse_mag_var2dd(self):
        if self.mag_var_src == '':  # No value - assume that bearings are true!
            self.is_valid = True
            self.mag_var_dd = 0.0
        else:
            mag_var_norm = self.normalize_src_input(self.mag_var_src)
            mag_var_valid, magvar_dd = self.check_angle_dd(mag_var_norm)
            if mag_var_valid is False:
                self.is_valid = False
                self.err_msg = 'Magnetic variation error!\n'
            else:
                self.mag_var_dd = magvar_dd

    @property
    def mag_var_dd(self):
        return self._mag_var_dd

    @mag_var_dd.setter
    def mag_var_dd(self, value):
        self._mag_var_dd = value


class Bearing(AngleBase):
    def __init__(self, brng_src, checked_value='Bearing'):
        AngleBase.__init__(self)
        self.brng_src = brng_src
        self._brng_dd = None
        self.checked_value = checked_value
        self.parse_brng2dd()

    def parse_brng2dd(self):
        """ Parse source value to convert it into decimal degrees value"""
        if self.brng_src == '':  # No value
            self.is_valid = False
            self.err_msg = 'Enter ' + self.checked_value + '\n'
        else:
            brng_norm = self.normalize_src_input(self.brng_src)
            # Check if angle is in DD format without hemisphere prefix, suffix
            self.is_valid, brng_dd = self.check_angle_dd(brng_norm)
            if brng_dd is None:
                self.is_valid, brng_dd = self.parse_compacted_formats(coord_lat_comp_regex, brng_norm)
                # Check separated format
                if brng_dd is None:
                    self.is_valid, brng_dd = self.parse_separated_formats(brng_norm)

            if brng_dd is not None:
                # Check baring range
                self.is_valid, self.brng_dd = self.check_angle_range(brng_dd, 0, 360)

            if self.is_valid is False:
                self.err_msg = self.checked_value + ' error\n'

    def calc_tbrng(self, mag_var_dd):
        """ Calculates true bearing.
        :param: dd_mag_var: float, magnetic variation value
        """
        if mag_var_dd == 0:
            return self.brng_dd
        else:
            tbrng = self.brng_dd + mag_var_dd
            if tbrng > 360:
                tbrng -= 360
            elif tbrng < 360:
                tbrng += 360
            return tbrng

    @property
    def brng_dd(self):
        return self._brng_dd

    @brng_dd.setter
    def brng_dd(self, value):
        self._brng_dd = value


""" Local coordinates """


class AviationBasePoint(AviationBaseClass):
    def __init__(self, ident='', lat_src='', lon_src='', mag_var_src=''):
        AviationBaseClass.__init__(self)
        self.ident = ident
        self.coordinates = CoordinatesPair(lat_src, lon_src)
        self.mag_var = MagVar(mag_var_src)
        self.check_base_point()

    def check_base_point(self):
        if self.coordinates.is_valid is False:
            self.is_valid = False
            self.err_msg += self.coordinates.err_msg

        if self.mag_var.is_valid is False:
            self.is_valid = False
            self.err_msg += self.mag_var.err_msg


class AviationCalculatedPoint(AngleBase):
    def __init__(self, ref_point: AviationBasePoint):
        AngleBase.__init__(self)
        self.ref_point = ref_point
        self._cp_lat_dd = None
        self._cp_lon_dd = None
        self._cp_definition = ''

    def polar_coordinates2latlon(self, angular_coord: Bearing, radial_coord: Distance):
        true_or_mag = 'TRUE'
        if self.ref_point.mag_var.mag_var_src != '':
            true_or_mag = 'MAG'

        tbrng = angular_coord.calc_tbrng(self.ref_point.mag_var.mag_var_dd)
        distance_m = radial_coord.get_meters()
        self.cp_lat_dd, self.cp_lon_dd = vincenty_direct_solution(self.ref_point.coordinates.lat_dd,
                                                                  self.ref_point.coordinates.lon_dd,
                                                                  tbrng,
                                                                  distance_m,
                                                                  WGS84_A, WGS84_B, WGS84_F)
        self.cp_definition = 'Ref. ident: {}; Azimuth: {} {}; Distance: {} {}'.format(self.ref_point.ident,
                                                                                      angular_coord.brng_src,
                                                                                      true_or_mag,
                                                                                      radial_coord.dist_src,
                                                                                      radial_coord.dist_uom)

    def offset_coordinates2latlon(self, brng: Bearing, dist: Distance, offset_side, offset_dist: Distance):
        tbrng = brng.calc_tbrng(self.ref_point.mag_var.mag_var_dd)
        dist_m = dist.get_meters()
        offset_dist_m = offset_dist.get_meters()

        self.cp_lat_dd, self.cp_lon_dd = dist_azm_orth_offset2latlon(self.ref_point.coordinates.lat_dd,
                                                                     self.ref_point.coordinates.lon_dd,
                                                                     tbrng,
                                                                     dist_m,
                                                                     offset_dist_m,
                                                                     offset_side)
        self.cp_definition = 'Ref. ident: {}; Distance: {} {}; Offset side: {} {} {}'.format(self.ref_point.ident,
                                                                                             dist.dist_src,
                                                                                             dist.dist_uom,
                                                                                             offset_side,
                                                                                             offset_dist.dist_src,
                                                                                             offset_dist.dist_uom)

    @staticmethod
    def calc_reverse_axis_brng(axis_brng):
        reverse_axis_brng = axis_brng - 180
        if reverse_axis_brng < 0:
            reverse_axis_brng += 360
        elif reverse_axis_brng > 360:
            reverse_axis_brng -= 360
        return reverse_axis_brng

    def local_cartesian_coordinates2latlon(self, x_axis_brng: Bearing, y_axis_orientation, x: Distance, y: Distance):
        x_axis_tbrng = x_axis_brng.calc_tbrng(self.ref_point.mag_var.mag_var_dd)
        if y_axis_orientation == 'LEFT':
            y_axis_tbrng = x_axis_tbrng - 90
        elif y_axis_orientation == 'RIGHT':
            y_axis_tbrng = x_axis_tbrng + 90

        if y_axis_tbrng > 360:
            y_axis_tbrng -= 360
        elif y_axis_tbrng < 0:
            y_axis_tbrng += 360

        reverse_x_axis_brng = self.calc_reverse_axis_brng(x_axis_tbrng)
        reverse_y_axis_brng = self.calc_reverse_axis_brng(y_axis_tbrng)

        x_m = x.get_meters()
        y_m = y.get_meters()

        if x_m >= 0:  # Normal x axis azimuth
            if y_m >= 0:  # Normal y axis azimuth
                self.cp_lat_dd, self.cp_lon_dd = dist_azm_orth_offset2latlon2(self.ref_point.coordinates.lat_dd,
                                                                              self.ref_point.coordinates.lon_dd,
                                                                              x_axis_tbrng,
                                                                              y_axis_tbrng,
                                                                              x_m, y_m)
            elif y_m < 0:
                self.cp_lat_dd, self.cp_lon_dd = dist_azm_orth_offset2latlon2(self.ref_point.coordinates.lat_dd,
                                                                              self.ref_point.coordinates.lon_dd,
                                                                              x_axis_tbrng,
                                                                              reverse_y_axis_brng,
                                                                              x_m, math.fabs(y_m))

        elif x_m < 0:
            if y_m >= 0:  # Normal y axis azimuth
                self.cp_lat_dd, self.cp_lon_dd = dist_azm_orth_offset2latlon2(self.ref_point.coordinates.lat_dd,
                                                                              self.ref_point.coordinates.lon_dd,
                                                                              reverse_x_axis_brng,
                                                                              y_axis_tbrng,
                                                                              math.fabs(x_m), y_m)
            elif y_m < 0:
                self.cp_lat_dd, self.cp_lon_dd = dist_azm_orth_offset2latlon2(self.ref_point.coordinates.lat_dd,
                                                                              self.ref_point.coordinates.lon_dd,
                                                                              reverse_x_axis_brng,
                                                                              reverse_y_axis_brng,
                                                                              math.fabs(x_m), math.fabs(y_m))

        self.cp_definition = 'Ref. ident: {}; X axis brng: {}; X: {} {}; Y: {} {}'.format(self.ref_point.ident,
                                                                                          x_axis_brng.brng_src,
                                                                                          x.dist_src, x.dist_uom,
                                                                                          y.dist_src, y.dist_uom)

    def get_calc_point_dms(self):
        """ Gets DMS format of calculated point
        :return: tuple: two elements of tuple, str """
        lat_dms = self.dd2dms(self.cp_lat_dd)
        if self.cp_lat_dd < 0:
            lat_suffix = 'S'
        else:
            lat_suffix = 'N'
        cp_lat_dms = lat_dms + lat_suffix

        lon_dms = self.dd2dms(self.cp_lon_dd)
        if self.cp_lon_dd < 0:
            lon_suffix = 'W'
        else:
            lon_suffix = 'E'

        cp_lon_dms = lon_dms + lon_suffix

        return cp_lat_dms, cp_lon_dms

    @property
    def cp_lat_dd(self):
        return self._cp_lat_dd

    @cp_lat_dd.setter
    def cp_lat_dd(self, value):
        self._cp_lat_dd = value

    @property
    def cp_lon_dd(self):
        return self._cp_lon_dd

    @cp_lon_dd.setter
    def cp_lon_dd(self, value):
        self._cp_lon_dd = value

    @property
    def cp_definition(self):
        return self._cp_definition

    @cp_definition.setter
    def cp_definition(self, value):
        self._cp_definition = value


VALID = 'VALID'
NOT_VALID = 'NOT_VALID'

def check_distance2(d):
    """ Distance validation. Uses float() function to check if parameters is a number
    :param d: string, distance to validate
    :return is_valid: True if distance is valid,
                     constant NOT_VALID if distance is not valid (e.g distance is less than 0)
    """
    try:
        dist = float(d)
        if dist < 0:  # Check if is less than 0
            dist = NOT_VALID
    except ValueError:
        dist = NOT_VALID
    return dist


def check_azm_dist(azm, dist):
    """ Checks if azimuth and distance in CSV file are correct.
    :param azm: azimuth to check
    :param dist: distance to check
    :return: is_valid: bool, True if input is valid, False otherwise
    :return: err_msg: str, string with error message
    """
    is_valid = True
    err_msg = ''
    a = Bearing(azm)

    if a.is_valid is False:
        is_valid = False
        err_msg += '*Azimuth value error*'
    if check_distance2(dist) == NOT_VALID:
        is_valid = False
        err_msg += '*Distance value error*'
    return is_valid, err_msg
