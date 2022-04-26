import re


AT_LAT = 'LAT'
AT_LON = 'LON'
AT_MAGVAR = 'MAGVAR'
AT_BEARING = 'BEARING'

COORDINATES = {
    AT_LON: {
        'DMSH_COMPACTED': re.compile(r'''(?P<deg>^180|^1[0-7]\d|^0\d{2})
                                         (?P<min>[0-5]\d)
                                         (?P<sec>[0-5]\d\.\d+|[0-5]\d)
                                         (?P<hem>[EW]$)
                                      ''', re.VERBOSE),
        'HDMS_COMPACTED': re.compile(r'''(?P<hem>^[EW])
                                         (?P<deg>180|1[0-7]\d|0\d{2})
                                         (?P<min>[0-5]\d)
                                         (?P<sec>[0-5]\d\.\d+$|[0-5]\d$)          
                                      ''', re.VERBOSE),
    },
    AT_LAT: {
        'DMSH_COMPACTED': re.compile(r'''(?P<deg>^90|^[0-8]\d)
                                         (?P<min>[0-5]\d)
                                         (?P<sec>[0-5]\d\.\d+|[0-5]\d)
                                         (?P<hem>[NS]$)
                                    ''', re.VERBOSE),
        'HDMS_COMPACTED': re.compile(r'''(?P<hem>^[NS])
                                         (?P<deg>90|[0-8]\d)
                                         (?P<min>[0-5]\d)
                                         (?P<sec>[0-5]\d\.\d+$|[0-5]\d$)
                                     ''', re.VERBOSE),
    }
}


class Coordinate:

    def __init__(self, ang_src, ang_type):
        self._ang_src = ang_src
        self._ang_type = ang_type
        self._ang_dd = None

    def _is_within_range(self):
        if self._ang_type == AT_LON:
            return -180 <= self._ang_dd <= 180
        elif self._ang_dd == AT_LAT:
            return -90 <= self._ang_dd <= 90

    @staticmethod
    def _convert_dmsh_parts_to_dd(d, m, s, h):
        """
        :param d: int, degrees
        :param m: int, minutes
        :param s: float, seconds
        :param h: str, hemisphere designator
        :return: float: decimal degrees
        """
        if (0 <= m < 60) and (0 <= s < 60):
            dd = d + m / 60 + s / 3600
            if h in ['W', 'S']:
                return -dd
            elif h in ['N', 'E']:
                return dd

    def convert_to_dd(self):
        norm_angle = self._ang_src
        for rgx in COORDINATES[self._ang_type].values():
            if rgx.match(norm_angle):
                groups = rgx.search(norm_angle)
                d = int(groups.group('deg'))
                m = int(groups.group('min'))
                s = float(groups.group('sec'))
                h = groups.group('hem')
                # TODO: move to Angle Base
                return Coordinate._convert_dmsh_parts_to_dd(d, m, s, h)
        else:
            raise ValueError('Error')



# import abc
# import re
# A_LAT = 'C_LAT'
# A_LON = 'C_LON'
#
# # Coordinate supported formats:
# # Compacted HDMS, DMSH
# # Separated HD-M-S, D-M-SH, HD M S, D M SH
#
#
#
# COORDINATE_COMPACTED = {
#     A_LON: {
#         'DMSH_COMPACTED': re.compile(r'''(?P<deg>^180|^1[0-7]\d|^0\d{2})
#                                          (?P<min>[0-5]\d)
#                                          (?P<sec>[0-5]\d\.\d+|[0-5]\d)
#                                          (?P<hem>[EW]$)
#                                       ''', re.VERBOSE),
#         'HDMS_COMPACTED': re.compile(r'''(?P<hem>^[EW])
#                                          (?P<deg>180|1[0-7]\d|0\d{2})
#                                          (?P<min>[0-5]\d)
#                                          (?P<sec>[0-5]\d\.\d+$|[0-5]\d$)
#                                       ''', re.VERBOSE),
#     },
#     A_LAT: {
#         'DMSH_COMPACTED': re.compile(r'''(?P<deg>^90|^[0-8]\d)
#                                          (?P<min>[0-5]\d)
#                                          (?P<sec>[0-5]\d\.\d+|[0-5]\d)
#                                          (?P<hem>[NS]$)
#                                     ''', re.VERBOSE),
#         'HDMS_COMPACTED': re.compile(r'''(?P<hem>^[NS])
#                                          (?P<deg>90|[0-8]\d)
#                                          (?P<min>[0-5]\d)
#                                          (?P<sec>[0-5]\d\.\d+$|[0-5]\d$)
#                                      ''', re.VERBOSE),
#     }
# }
#
# from abc import ABC, abstractmethod
# from unittest.mock import Mock, MagicMock
#
#
# class AbstractAngle(ABC):
#
#     @abstractmethod
#     def __init__(self):
#         ABC.__init__(self)
#
#     # @property
#     # @abstractmethod
#     # def ang_dd(self):
#     #     pass
#
#     @abstractmethod
#     def _is_within_range(self, *args):
#         pass
#
#     @abstractmethod
#     def convert_to_dd(self):
#         pass
#
#
# class MagVar(AbstractAngle):
#
#     def __init__(self, ang_src):
#         AbstractAngle.__init__(self)
#         self._ang_src = ang_src
#         self._ang_dd = None
#         self.convert_to_dd()
#
#     def _is_within_range(self):
#         """
#
#         :return:
#         :rtype:
#         >>> MagVar = mock()
#         >>> MagVar._ang_dd = 100
#         >>> assert MagVar()._is_within_range() is False
#         >>> mock = MagicMock()
#         >>> mock.__str__.return_value = 'foobarbaz'
#         >>> assert MagVar(100)._is_within_range() is False
#         """
#         if self._ang_dd:
#             return -90 <= self._ang_dd <= 90
#
#     def convert_to_dd(self):
#         self._ang_dd = self._ang_src
#
#
#
# # class BaseAngle:
# #
# #     def __init__(self, ang_src):
# #         self._ang_src = ang_src
# #         self._ang_dd = None
# #
# #     @staticmethod
# #     def _check_angle_range(angle_dd, angle_type):
# #         """ Check if angle is within closed interval <min_value, max_value>
# #         :param angle_dd: float, angle value to check
# #         :param angle_type: str
# #         :return: bool, True if angle is within the range, False otherwise
# #         """
# #         if angle_type == A_LON:
# #             return -180 <= angle_dd <= 180
# #         elif angle_type == A_LAT:
# #             return -90 <= angle_dd <= 90
# #
# #
# # class BaseAngle2:
# #
# #     def __init__(self):
# #         pass
# #
# #     @staticmethod
# #     def _check_angle_range(angle_dd, angle_type):
# #         """ Check if angle is within closed interval <min_value, max_value>
# #         :param angle_dd: float, angle value to check
# #         :param angle_type: str
# #         :return: bool, True if angle is within the range, False otherwise
# #         >>> assert BaseAngle._check_angle_range(-180.00001, A_LON) is False
# #         >>> assert BaseAngle._check_angle_range(-180.0000, A_LON) is True
# #         >>> assert BaseAngle._check_angle_range(-179.999, A_LON) is True
# #         >>> assert BaseAngle._check_angle_range(179.999, A_LON) is True
# #         >>> assert BaseAngle._check_angle_range(180.000, A_LON) is True
# #         >>> assert BaseAngle._check_angle_range(180.0001, A_LON) is False
# #         >>> assert BaseAngle._check_angle_range(-90.00001, A_LAT) is False
# #         >>> assert BaseAngle._check_angle_range(-90.0000, A_LAT) is True
# #         >>> assert BaseAngle._check_angle_range(-89.999, A_LAT) is True
# #         >>> assert BaseAngle._check_angle_range(89.999, A_LAT) is True
# #         >>> assert BaseAngle._check_angle_range(90.000, A_LAT) is True
# #         >>> assert BaseAngle._check_angle_range(90.0001, A_LAT) is False
# #         """
# #         if angle_type == A_LON:
# #             return -180 <= angle_dd <= 180
# #         elif angle_type == A_LAT:
# #             return -90 <= angle_dd <= 90
# #
# #     @staticmethod
# #     def _convert_dmsh_parts_to_dd(d, m, s, h):
# #         """
# #         :param d: int, degrees
# #         :param m: int, minutes
# #         :param s: float, seconds
# #         :param h: str, hemisphere designator
# #         :return: float: decimal degrees
# #         >>> assert BaseAngle._convert_dmsh_parts_to_dd(100, 35, 44, 'N') == 100.59555555555555
# #         >>> assert BaseAngle._convert_dmsh_parts_to_dd(100, 35, 44, 'W') == -100.59555555555555
# #         >>> assert BaseAngle._convert_dmsh_parts_to_dd(100, 0, 0, 'A') is None
# #         >>> assert BaseAngle._convert_dmsh_parts_to_dd(100, -1, 0, 'E') is None
# #         >>> assert BaseAngle._convert_dmsh_parts_to_dd(100, 61, 0, 'E') is None
# #         >>> assert BaseAngle._convert_dmsh_parts_to_dd(100, 0, -1, 'E') is None
# #         >>> assert BaseAngle._convert_dmsh_parts_to_dd(100, 0, -0.001, 'E') is None
# #         >>> assert BaseAngle._convert_dmsh_parts_to_dd(100, 0, 60.001, 'E') is None
# #         """
# #         if (0 <= m < 60) and (0 <= s < 60):
# #             dd = d + m / 60 + s / 3600
# #             if h in ['W', 'S']:
# #                 return -dd
# #             elif h in ['N', 'E']:
# #                 return dd
# #
# #     @staticmethod
# #     def _compacted_format_to_dd(norm_angle, angl_type):
# #         """
# #         :param norm_angle: str
# #         :param angl_type: str
# #         >>> assert BaseAngle._compacted_format_to_dd('1800000.0E', A_LON) == 180.0
# #         >>> assert BaseAngle._compacted_format_to_dd('1800000.0W', A_LON) == -180.0
# #         """
# #         for pattern_name, pattern_rgx in COORDINATE_COMPACTED[angl_type].items():
# #             if pattern_rgx.match(norm_angle):
# #                 if 'H' in pattern_name:
# #                     groups = pattern_rgx.search(norm_angle)
# #                     d = int(groups.group('deg'))
# #                     m = int(groups.group('min'))
# #                     s = float(groups.group('sec'))
# #                     h = groups.group('hem')
# #                     return BaseAngle._convert_dmsh_parts_to_dd(d, m, s, h)
# #
# #
# # class Coordinate(BaseAngle):
# #
# #     def __init__(self, coord_src, coord_type, label):
# #         """
# #         >>> valid_longitudes = [('1800000E', 180.0),
# #         ...                     ('1003000.0E', 100.5),
# #         ...                     ('0103000.0E', 10.5),
# #         ...                     ('0013000.0E', 1.5),
# #         ...                     ('0003000.0E', 0.5),
# #         ...                     ('1234601.445E', 123.76706805555555),
# #         ...                     ('1800000W', -180.0),
# #         ...                     ('1003000.0W', -100.5),
# #         ...                     ('0103000.0W', -10.5),
# #         ...                     ('0013000.0W', -1.5),
# #         ...                     ('0003000.0W', -0.5),
# #         ...                     ('E1800000', 180.0),
# #         ...                     ('E1003000.0', 100.5),
# #         ...                     ('E0103000.0', 10.5),
# #         ...                     ('E0013000.0', 1.5),
# #         ...                     ('E0003000.0', 0.5),
# #         ...                     ('E1234601.445', 123.76706805555555),
# #         ...                     ('W1800000', -180.0),
# #         ...                     ('W1003000.0', -100.5),
# #         ...                     ('W0103000.0', -10.5),
# #         ...                     ('W0013000.0', -1.5),
# #         ...                     ('W0003000.0', -0.5),
# #         ...                     ('W1234601.445', -123.76706805555555)]
# #         >>> for lon_src, lon_dd in valid_longitudes:
# #         ...    c = Coordinate(lon_src, A_LON, 'Longitude')
# #         ...    assert c.coord_dd == lon_dd, 'actual: {}, expected: {}'.format(c.coord_dd, lon_dd)
# #         >>> invalid_longitudes = ['1800000.01E']
# #         >>> for lon_src in invalid_longitudes:
# #         ...   c = Coordinate(lon_src, A_LON, 'Longitude')
# #         Traceback (most recent call last):
# #         ValueError: Cannot convert 1800000.01E to DD.
# #         <BLANKLINE>
# #         """
# #         BaseAngle.__init__(self)
# #         self._coord_src = coord_src
# #         self._coord_type = coord_type
# #         self._label = label
# #         self._coord_dd = None
# #         self._convert_to_dd()
# #
# #     @property
# #     def coord_dd(self):
# #         return self._coord_dd
# #
# #     def _convert_to_dd(self):
# #         """
# #
# #         """
# #         coord_norm = self._coord_src
# #         dd = Coordinate._compacted_format_to_dd(coord_norm, self._coord_type)
# #         if dd is None:
# #             raise ValueError('Cannot convert {} to DD.\n'.format(self._coord_src))
# #
# #         if not BaseAngle._check_angle_range(dd, self._coord_type):
# #             raise ValueError('Cannot convert {} to DD.\n'.format(self._coord_src))
# #
# #         self._coord_dd = dd
# #
# #
# # # Coordinate('1800000.01E', A_LON, 'Longitude')
# #
# # class Point:
# #
# #     def __init__(self, label, lon_src, lat_src, mag_var=None):
# #         try:
# #             self._lon = Coordinate(lon_src, coord_type=A_LON, label=label)
# #         except ValueError as e:
# #             self._err_msg += e
# #
# #         try:
# #             self._lat = Coordinate(lat_src, coord_type=A_LON, label=label)
# #         except ValueError:
# #             self._err_msg += e
# #
# #         try:
# #             self._mag_var = MagVar(mag_var) if mag_var else
#
