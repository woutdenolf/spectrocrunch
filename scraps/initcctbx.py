# -*- coding: utf-8 -*-

"""Initialize cctbx library. To use it, add to following line to the top of your script:
execfile("initcctbx.py")
"""

CCTBCBUILD = "/usr/local/cctbx/build"

import os, sys

os.environ["LIBTBX_BUILD"] = CCTBCBUILD

FONTCONFIG_PATH = os.path.abspath(
    os.path.join(CCTBCBUILD, "..", "base", "etc", "fonts")
)
os.environ["FONTCONFIG_PATH"] = FONTCONFIG_PATH

os.environ["LD_LIBRARY_PATH"] = ":".join(
    [
        os.path.abspath(os.path.join(CCTBCBUILD, "lib")),
        os.path.abspath(os.path.join(CCTBCBUILD, "bin")),
        os.path.abspath(os.path.join(CCTBCBUILD, "..", "base", "lib")),
        os.environ["LD_LIBRARY_PATH"],
    ]
)
try:
    os.execl(sys.executable, "python", __file__, *sys.argv[1:])
    sys.exit()
except Exception as exc:
    print("Failed re-exec:", exc)
    sys.exit(1)

sys.path.append(os.path.abspath(os.path.join(CCTBCBUILD, "lib")))
sys.path.append(
    os.path.abspath(
        os.path.join(
            CCTBCBUILD, "..", "modules", "cctbx_project", "libtbx", "pythonpath"
        )
    )
)
sys.path.append(
    os.path.abspath(
        os.path.join(CCTBCBUILD, "..", "modules", "cctbx_project", "boost_adaptbx")
    )
)
sys.path.append(os.path.abspath(os.path.join(CCTBCBUILD, "..", "modules")))
sys.path.append(
    os.path.abspath(os.path.join(CCTBCBUILD, "..", "modules", "cctbx_project"))
)
