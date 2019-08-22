# -*- coding: utf-8 -*-

from __future__ import absolute_import
import warnings

try:
    import xraylib_np
except ImportError:
    xraylib_np = None

try:
    from xraylib import *
except ImportError:
    XRayInit = None
    warnings.warn("xraylib is not installed", ImportWarning)
else:
    XRayInit()
    SetErrorMessages(0)

    # Code <-> Name: one-to-one
    code_to_shell = {v: s.split(
        '_')[0] for s, v in globals().items() if s.endswith("_SHELL")}
    shell_to_code = {v: k for k, v in code_to_shell.items()}
    shell_min = min(code_to_shell)
    shell_max = max(code_to_shell)

    # Code <-> Name: one-to-many
    line_to_code = {s.split('_')[0]: v for s,
                    v in globals().items() if s.endswith("_LINE")}
    code_to_line = {code: [name for name, code2 in line_to_code.items(
    ) if code == code2] for code in set(line_to_code.values())}
    line_min = min(code_to_line)
    line_max = max(code_to_line)

    # Composites
    composites = {KA_LINE: [KL3_LINE, KL2_LINE, KL1_LINE],
                  KB_LINE: [KP5_LINE, KP4_LINE, KP3_LINE, KP2_LINE, KP1_LINE] +
                  [KO7_LINE, KO6_LINE, KO5_LINE, KO4_LINE, KO3_LINE, KO2_LINE, KO1_LINE] +
                  [KN7_LINE, KN6_LINE, KN5_LINE, KN4_LINE, KN3_LINE, KN2_LINE, KN1_LINE] +
                  [KM5_LINE, KM4_LINE, KM3_LINE,
                   KM2_LINE, KM1_LINE],
                  KP_LINE: [KP5_LINE, KP4_LINE, KP3_LINE, KP2_LINE, KP1_LINE],
                  KO_LINE: [KO7_LINE, KO6_LINE, KO5_LINE, KO4_LINE, KO3_LINE, KO2_LINE, KO1_LINE],
                  LA_LINE: [L3M5_LINE, L3M4_LINE],
                  LB_LINE: [L3O4_LINE, L3O3_LINE, L3N5_LINE, L3N1_LINE, L2M4_LINE, L1M3_LINE, L1M2_LINE],
                  L1N67_LINE: [L1N7_LINE, L1N6_LINE],
                  L1O45_LINE: [L1O5_LINE, L1O4_LINE],
                  L1P23_LINE: [L1P3_LINE, L1P2_LINE],
                  L2P23_LINE: [L2P3_LINE, L2P2_LINE],
                  L3O45_LINE: [L3O5_LINE, L3O4_LINE],
                  L3P23_LINE: [L3P3_LINE, L3P2_LINE],
                  L3P45_LINE: [L3P5_LINE, L3P4_LINE]
                  }

    rcomposites = {}
    for k, v in composites.items():
        for l in v:
            if l in rcomposites:
                rcomposites[l].append(k)
                rcomposites[l].sort()
            else:
                rcomposites[l] = [k]
