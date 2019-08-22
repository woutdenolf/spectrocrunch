# -*- coding: utf-8 -*-

from . import element


def elemental_crosssections(dictin, dictout=None, w=1):
    if dictout is None:
        dictout = {}
    for k, v in dictin.items():
        if isinstance(v["cs"], dict):
            elemental_crosssections(v["cs"], w=v["w"], dictout=dictout)
        else:
            if k in dictout:
                dictout[k] += w*v["w"]*v["cs"]
            else:
                dictout[k] = w*v["w"]*v["cs"]
    return dictout
