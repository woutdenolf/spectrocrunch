import collections
from array import array
import numpy as np
from uncertainties import ufloat
from ...patch.pint import ureg


def gendata():
    data = {}
    data["type"] = type
    data["none"] = None
    data["bool"] = True
    data["int"] = int(1)
    data["float"] = float(1)
    data["npint"] = np.int16(1)
    data["npfloat"] = np.float16(1)
    data["nan"] = float("nan")
    data["inf"] = float("inf")
    data["str"] = "abc"
    data["unicode"] = "\u0101bc"
    data["bytes"] = b"\xe4bc"
    data["npunicode"] = np.str_("\u0101bc")
    data["npbytes"] = np.bytes_(b"\xe4bc")
    data["npstr"] = np.str_("abc")
    data["unicode_list"] = ["\u0101bc", "\u0101bc"]
    data["bytes_list"] = [b"\xe4bc", b"\xe4bc"]
    data["str_list"] = ["abc", "abc"]
    data["npunicode_array"] = np.array(["\u0101bc", "\u0101bc"])
    data["npbytes_array"] = np.array([b"\xe4bc", b"\xe4bc"])
    data["npstr_array"] = np.array(["abc", "abc"])
    data["bytearray"] = bytearray([228, 50, 51])
    data["array"] = array("i", [1000, 2000])
    tmp = 1, 2
    data["tuple"] = tuple(tmp)
    data["list"] = list(tmp)
    try:
        data["xrange"] = xrange(5)
    except NameError:
        data["xrange"] = range(5)
    data["frozenset"] = frozenset(tmp)
    data["set"] = set(tmp)
    data["ndarray"] = np.zeros(5)
    data["ndarray0"] = np.array(5)

    class Dummy:
        pass

    data["nparray"] = np.array([Dummy()] * 5)
    data["nparray0"] = np.array(Dummy())
    data["qarray"] = ureg.Quantity(np.zeros(5), ureg.meter)
    data["qarray0"] = ureg.Quantity(np.array(5), ureg.meter)
    data["qscalar"] = ureg.Quantity(np.float32(5), ureg.meter)
    tmp = [("yellow", 1), ("blue", 2)]
    data["dict"] = dict(tmp)
    data["odict"] = collections.OrderedDict(tmp)
    data["ddict"] = collections.defaultdict(tuple, tmp)
    data["uscalar"] = ufloat(1, 2)
    data["uarray"] = np.array([ufloat(1, 2), ufloat(3, 4)])
    data["uarray0"] = np.array(ufloat(1, 2))
    return data
