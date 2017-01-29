# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# Imports
import sys
import subprocess
import os
from setuptools import setup, Command, find_packages
from setuptools.command.install import install
import _version

# Disable hardlinks when not working
if hasattr(os, 'link'):
    tempfile = __file__ + '.tmp'
    try:
        os.link(__file__, tempfile)
    except OSError as e:
        if e.errno == 1:  # Operation not permitted
            del os.link
        else:
            raise
    finally:
        if os.path.exists(tempfile):
            os.remove(tempfile)

# Get setup information
def get_version():
    return _version.strictversion

def get_devstatus():
    # The development status is derived from the SpectroCrunch release level
    mapping = {"dev":2,"alpha":3,"beta":4,"rc":5,"final":6}
    cycle = {1:"Planning",2:"Pre-Alpha",3:"Alpha",4:"Beta",5:"Production/Stable",6:"Mature",7:"Inactive"}

    status = mapping[_version.version_info.releaselevel]
    
    return "Development Status :: %d - %s"%(status,cycle[status])

def get_readme():
    dirname = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(dirname, "README.rst"), "r") as fp:
        long_description = fp.read()
    return long_description

# Command classes
cmdclass = {}

class TestAllPackages(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):    
        errno = subprocess.call([sys.executable,'-m','spectrocrunch.tests.test_all'])
        if errno != 0:
            print("Tests did not pass !!!")
            raise SystemExit(errno)
        else:
            print("All Tests passed.")
cmdclass['test'] = TestAllPackages

class VersionOfAllPackages(Command):
    user_options = []
    
    def initialize_options(self):
        pass
    
    def finalize_options(self):
        pass
    
    def run(self):
        print("This version of SpectroCrunch is", _version.version)
cmdclass['version'] = VersionOfAllPackages

class InstallWithVersion(install):
    def run(self):
        install.run(self)
cmdclass['install'] = InstallWithVersion

# Trove classifiers
classifiers = [get_devstatus(),
               "Environment :: Console",
               #"Environment :: MacOS X",
               #"Environment :: Win32 (MS Windows)",
               #"Environment :: X11 Applications :: Qt",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Natural Language :: English",
               "Operating System :: Microsoft :: Windows",
               "Operating System :: POSIX :: Linux",
               #"Operating System :: MacOS :: MacOS X",
               "Programming Language :: Python :: 2.7",
               "Programming Language :: Python :: 3.4",
               "Programming Language :: Python :: 3.5",
               #"Topic :: Documentation :: Sphinx",
               "Topic :: Scientific/Engineering :: Physics",
               "Topic :: Software Development :: Libraries :: Python Modules",
               ]

# Needed for using Spectrocrunch
install_requires = ["numpy", "scipy", "h5py", "fabio", "silx", "pyparsing", "PyMca5", "shapely", "matplotlib"]
extras_require = {\
    "physics":["xraylib", "cctbx", "fdmnes"],\
    "elastix":["SimpleITK"]\
    }

# Needed for running the setup script
setup_requires = ["testfixtures"]

setup(name='SpectroCrunch',
      version=get_version(),
      url="https://github.com/woutdenolf/spectrocrunch",
      author="Wout De Nolf",
      author_email="woutdenolf@users.sf.net",
      classifiers = classifiers,
      description="Spectroscopy data crunching",
      long_description=get_readme(),
      packages=find_packages(),
      install_requires=install_requires,
      extras_require=extras_require,
      setup_requires=setup_requires,
      include_package_data=True,
      license="MIT",
      cmdclass=cmdclass
      )
