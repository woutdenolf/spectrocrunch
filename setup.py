# -*- coding: utf-8 -*-

import _version
from setuptools.command.build_py import build_py
from setuptools.command.install import install
from setuptools import Command
from setuptools import find_packages
from setuptools import setup
import fnmatch
import glob
import shutil
import os
import subprocess
import sys

try:
    import sphinx
    from sphinx.setup_command import BuildDoc
except ImportError:
    sphinx = None


PROJECT = "spectrocrunch"

########################################
## Disable hardlinks when not working ##
########################################
if hasattr(os, "link"):
    tempfile = __file__ + ".tmp"
    try:
        os.link(__file__, tempfile)
    except OSError as e:
        del os.link
    finally:
        if os.path.exists(tempfile):
            os.remove(tempfile)


###########################
## Get setup information ##
###########################
def get_version():
    return _version.strictversion


def get_devstatus():
    # The development status is derived from the SpectroCrunch release level
    mapping = {"dev": 2, "alpha": 3, "beta": 4, "rc": 5, "final": 6}
    cycle = {
        1: "Planning",
        2: "Pre-Alpha",
        3: "Alpha",
        4: "Beta",
        5: "Production/Stable",
        6: "Mature",
        7: "Inactive",
    }

    status = mapping[_version.version_info.releaselevel]

    return "Development Status :: %d - %s" % (status, cycle[status])


def get_readme():
    dirname = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(dirname, "README.rst"), "r") as fp:
        long_description = fp.read()
    return long_description


#####################
## Command classes ##
#####################
cmdclass = {}


class DisabledCommand(Command):
    user_options = []

    _MSG = "Command is disabled."

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        raise RuntimeError(self._MSG)


#######################
## "version" command ##
#######################
class VersionOfAllPackages(Command):
    description = "Get project version"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        print("This version of {} is {}".format(PROJECT, _version.version))


cmdclass["version"] = VersionOfAllPackages


########################
## "build_py" command ##
########################
class BuildWithVersion(build_py):
    """
    Enhanced build_py which copies version.py to <PROJECT>._version.py
    """

    description = "build with version info"

    def find_package_modules(self, package, package_dir):
        modules = build_py.find_package_modules(self, package, package_dir)
        if "." not in package:
            modules.append((package, "_version", "_version.py"))
        return modules


cmdclass["build_py"] = BuildWithVersion


#########################
## "build_doc" command ##
#########################
if sphinx is not None:

    class BuildDocCommand(BuildDoc):
        description = "Build documentation from source"

        def run(self):
            # Make sure the python path is pointing to the newly built
            # code so that the documentation is built on this and not a
            # previously installed version
            build = self.get_finalized_command("build")
            sys.path.insert(0, os.path.abspath(build.build_lib))

            for builder in ["html", "latex"]:
                self.builder = builder
                self.builder_target_dir = os.path.join(self.build_dir, builder)
                self.mkpath(self.builder_target_dir)
                BuildDoc.run(self)

            sys.path.pop(0)

else:

    class BuildDocCommand(DisabledCommand):
        _MSG = "Sphinx is required to build or test the documentation."


cmdclass["build_doc"] = BuildDocCommand


#####################
## "clean" command ##
#####################
class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""

    description = "Clean build and compiled files"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        shutil.rmtree("./build", True)
        shutil.rmtree("./dist", True)

        patterns = ["*.egg-info"]
        for pattern in patterns:
            for dirname in glob.glob(pattern):
                shutil.rmtree(dirname, True)

        patterns = ["*.pyc"]
        for root, dirnames, filenames in os.walk(PROJECT):
            for pattern in patterns:
                for filename in fnmatch.filter(filenames, pattern):
                    os.remove(os.path.join(root, filename))


cmdclass["clean"] = CleanCommand


#####################
## "name" command ##
#####################
class NameCommand(Command):
    """Print project name."""

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        print(PROJECT)


cmdclass["name"] = NameCommand


#######################
## Trove classifiers ##
#######################
classifiers = [
    get_devstatus(),
    "Environment :: Console",
    ## 'Environment :: MacOS X',
    ## 'Environment :: Win32 (MS Windows)',
    ## 'Environment :: X11 Applications :: Qt',
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Natural Language :: English",
    ## 'Operating System :: Microsoft :: Windows',
    "Operating System :: POSIX :: Linux",
    ## 'Operating System :: MacOS :: MacOS X',
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Topic :: Documentation :: Sphinx",
    "Topic :: Scientific/Engineering :: Physics",
    "Topic :: Software Development :: Libraries :: Python Modules",
]


##################
## Requirements ##
##################
install_requires = [
    "setuptools",
    "numpy",
    "future",
    "scipy",
    "h5py",
    "fabio",
    "silx",
    "pyparsing",
    "shapely",
    "matplotlib",
    "uncertainties",
    "pint<0.20",
    "pandas",
    "scikit-image",
    "xlsxwriter",
    "xlrd",
    "openpyxl",
    "python-dateutil",
    "jsonpickle",
    "testfixtures",
    "future",
    "cvxopt",
    "pymca5",
]
extras_require = {
    "physics": ["xraylib", "cctbx", "fdmnes", "PyTMM"],
    "elastix": ["SimpleITK"],
}
setup_requires = ["setuptools", "testfixtures"]


###################
## Package setup ##
###################
setup(
    name=PROJECT,
    version=get_version(),
    url="https://github.com/woutdenolf/spectrocrunch",
    author="Wout De Nolf",
    author_email="woutdenolf@users.sf.net",
    classifiers=classifiers,
    description="Spectroscopic imaging library (XRF/XAS)",
    long_description=get_readme(),
    install_requires=install_requires,
    extras_require=extras_require,
    setup_requires=setup_requires,
    packages=find_packages(),
    package_data={"spectrocrunch.resources": ["*/*.*"]},
    license="MIT",
    cmdclass=cmdclass,
    test_suite="{}.tests.test_all.test_suite".format(PROJECT),
)
