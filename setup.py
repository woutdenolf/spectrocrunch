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
"""
Package installation file to be used as:

- python setup.py version
- python setup.py test
- python setup.py build
- python setup.py install
- python setup.py sdist
- python setup.py bdist

"""

# Imports
from __future__ import print_function
import sys
import os
import subprocess
from setuptools import setup, Command, find_packages

# Get setup information
def get_version():
    import imp
    _version = imp.load_source('_version', os.path.join(os.path.dirname(os.path.abspath(__file__)), "_version.py"))    
    return _version.strictversion

def get_devstatus():
    # The development status is derived from the SpectroCrunch release level
    mapping = {"dev":2,"alpha":3,"beta":4,"rc":5,"final":6}
    cycle = {1:"Planning",2:"Pre-Alpha",3:"Alpha",4:"Beta",5:"Production/Stable",6:"Mature",7:"Inactive"}

    import imp
    _version = imp.load_source('_version', os.path.join(os.path.dirname(os.path.abspath(__file__)), "_version.py"))   

    status = mapping[_version.version_info.releaselevel]
    
    return "Development Status :: %d - %s"%(status,cycle[status])


def get_readme():
    dirname = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(dirname, "README.rst"), "r") as fp:
        long_description = fp.read()
    return long_description

# Command class for testing
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
            # raise SystemExit(errno)
        else:
            print("All Tests passed.")

class VersionOfAllPackages(Command):
    user_options = []
    
    def initialize_options(self):
        pass
    
    def finalize_options(self):
        pass
    
    def run(self):
        print("This version of SpectroCrunch is", get_version())


# Setup
cmdclass = {'test':TestAllPackages,'version':VersionOfAllPackages}

# Trove classifiers
classifiers = [get_devstatus(),
               "Environment :: Console",
               "Environment :: MacOS X",
               "Environment :: Win32 (MS Windows)",
               "Environment :: X11 Applications :: Qt",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Natural Language :: English",
               "Operating System :: Microsoft :: Windows",
               "Operating System :: POSIX",
               "Operating System :: MacOS :: MacOS X",
               "Programming Language :: Python :: 2",
               "Programming Language :: Python :: 3",
               "Topic :: Documentation :: Sphinx",
               "Topic :: Scientific/Engineering :: Physics",
               "Topic :: Software Development :: Libraries :: Python Modules",
               ]

# Needed for using Spectrocrunch
install_requires = ["numpy", "scipy", "h5py", "fabio", "sift_pyocl", "pyparsing", "PyMca5"]
# Non pip: "xraylib", "cctbx", "fdmnes", "SimpleITK"

# Needed for running the setup script
setup_requires = ["testfixtures", "setuptools_scm"] 

setup(name='SpectroCrunch',
      version=get_version(),
      #url="https://github.com/spectrocrunch",
      author="Wout De Nolf",
      author_email="woutdenolf@users.sf.net",
      classifiers = classifiers,
      description="Spectroscopy data crunching",
      long_description=get_readme(),
      packages=find_packages(),
      install_requires=install_requires,
      setup_requires=setup_requires,
      include_package_data=True,
      use_scm_version=True,
      license="MIT",
      cmdclass=cmdclass
      )
