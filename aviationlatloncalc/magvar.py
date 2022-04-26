"""
Supported formats:
 - DD: decimal degrees
 - DDH: decimal degrees with hemisphere suffix , example: 11.4E
 - HDD: decimal degrees with hemisphere prefix, example: E11.4
"""
from angle import BaseAngle
import re


MAGVAR_REGEXES = {
    'HDD': re.compile(r'''(?P<hem>^[EW])
                          (?P<dd>\d{1,2}\.\d+$|\d{1,2}$)
                      ''', re.VERBOSE),
    'DDH': re.compile(r'''(?P<dd>^\d{1,2}\.\d+|^\d{1,2})
                          (?P<hem>[EW]$)
                      ''', re.VERBOSE)
}


class MagVar(BaseAngle):

    """
    Attributes
    ----------
    _magvar_src: str
        source value of magnetic variation
    _magvar_dd: float
        magnetic variation converted to decimal degrees
    """

    def __init__(self, magvar_src=None):
        """
        >>> valid_mag_vars = [(11, 11),
        ...                   (-11, -11),
        ...                   (11.4, 11.4),
        ...                   (-11.4, -11.4),
        ...                   ('E1', 1.0),
        ...                   ('E1.123', 1.123),
        ...                   ('E11', 11.0),
        ...                   ('E11.4', 11.4),
        ...                   ('1E', 1.0),
        ...                   ('1.123E', 1.123),
        ...                   ('11E', 11.0),
        ...                   ('11.4E', 11.4),
        ...                   ('W1', -1.0),
        ...                   ('W1.123', -1.123),
        ...                   ('W11', -11.0),
        ...                   ('W11.4', -11.4),
        ...                   ('1W', -1.0),
        ...                   ('1.123W', -1.123),
        ...                   ('11W', -11.0),
        ...                   ('11.4W', -11.4)]
        >>> for mag_var_src, mag_var_dd in valid_mag_vars:
        ...    m = MagVar(mag_var_src)
        ...    assert m.magvar_dd == mag_var_dd, 'actual: {}, expected: {}'.format(m.magvar_dd, mag_var_dd)

        >>> invalid_mag_vars = [-100, 100, -100.1, 100.1, 'E100', 'W100', '100E', '100W', 'A10', '10A']

        """
        BaseAngle.__init__(self)
        self._magvar_src = magvar_src
        self._magvar_dd = None
        self._convert_to_dd()
        self._err_msg = ''
        if self._magvar_dd is None:
            self._err_msg = 'Mag var error: value {} is invalid or in not supported format.'.format(self._magvar_src)

    @property
    def magvar_dd(self):
        return self._magvar_dd

    @property
    def err_msg(self):
        return self._err_msg

    def __str__(self):
        return 'Mag var: {}'.format(self._magvar_src)

    def _normalize(self):
        return self._magvar_src

    def _convert_to_dd(self):
        # Check decimal degrees format
        if isinstance(self._magvar_src, (int, float)):
            if -90 <= self._magvar_src <= 90:
                self._magvar_dd = self._magvar_src
                return
        else:  # Check HDD, DDH formats
            norm_src = self._normalize()
            for rgx in MAGVAR_REGEXES.values():
                if rgx.match(norm_src):
                    groups = rgx.search(norm_src)
                    dd = float(groups.group('dd'))
                    h = groups.group('hem')

                    self._magvar_dd = -dd if h == 'W' else dd

# w Point -> laaie yjaktow i na koniec jeden wyjatek


a = MagVar(-101)
print(a.err_msg)
# MagVar(100)
# MagVar(-100.1)
# MagVar(100.1)