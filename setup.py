## -*- coding: utf-8 -*-
#
##   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
#
##   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
#
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
#
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
## THE SOFTWARE.

PROJECT = 'spectrocrunch'

#############
## Imports ##
#############
import sys
import subprocess
import os

from setuptools import setup
from setuptools import find_packages
from setuptools import Command
from setuptools.command.install import install
from setuptools.command.build_py import build_py

import _version

try:
    import sphinx
    from sphinx.setup_command import BuildDoc
except ImportError:
    sphinx = None


########################################
## Disable hardlinks when not working ## 
########################################
if hasattr(os, 'link'):
    tempfile = __file__ + '.tmp'
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
    ## # The development status is derived from the SpectroCrunch release level
    mapping = {"dev":2,"alpha":3,"beta":4,"rc":5,"final":6}
    cycle = {1:"Planning",2:"Pre-Alpha",3:"Alpha",4:"Beta",5:"Production/Stable",6:"Mature",7:"Inactive"}

    status = mapping[_version.version_info.releaselevel]
    
    return "Development Status :: %d - %s"%(status,cycle[status])

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


####################
## "test" command ## 
####################
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


#######################
## "version" command ## 
#######################
class VersionOfAllPackages(Command):
    user_options = []
    
    def initialize_options(self):
        pass
    
    def finalize_options(self):
        pass
    
    def run(self):
        print("This version of {} is {}".format(PROJECT,_version.version))
cmdclass['version'] = VersionOfAllPackages


########################
## "build_py" command ## 
########################
class BuildWithVersion(build_py):
    """
    Enhanced build_py which copies version.py to <PROJECT>._version.py
    """
    def find_package_modules(self, package, package_dir):
        modules = build_py.find_package_modules(self, package, package_dir)
        if "." not in package:
            modules.append((package, '_version', '_version.py'))
        return modules

cmdclass['build_py'] = BuildWithVersion


#########################
## "build_doc" command ## 
#########################
if sphinx is not None:
    class BuildDocCommand(BuildDoc):

        def run(self):
            ## # make sure the python path is pointing to the newly built
            ## # code so that the documentation is built on this and not a
            ## # previously installed version

            build = self.get_finalized_command('build')
            sys.path.insert(0, os.path.abspath(build.build_lib))

##             ## # Copy .ui files to the path:
##             dst = os.path.join(
##                 os.path.abspath(build.build_lib), "silx", "gui")
##             if not os.path.isdir(dst):
##                 os.makedirs(dst)
##             for i in os.listdir("gui"):
##                 if i.endswith(".ui"):
##                     src = os.path.join("gui", i)
##                     idst = os.path.join(dst, i)
##                     if not os.path.exists(idst):
##                         shutil.copy(src, idst)

            ## # Build the Users Guide in HTML and TeX format
            for builder in ['html', 'latex']:
                self.builder = builder
                self.builder_target_dir = os.path.join(self.build_dir, builder)
                self.mkpath(self.builder_target_dir)
                BuildDoc.run(self)
            sys.path.pop(0)
else:
    class BuildDocCommand(DisabledCommand):
        _MSG = 'Sphinx is required to build or test the documentation.'

cmdclass['build_doc'] = BuildDocCommand

#####################
## "clean" command ## 
#####################
class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')

cmdclass['clean'] = CleanCommand

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
        print PROJECT

cmdclass['name'] = NameCommand


#######################
## Trove classifiers ## 
#######################
classifiers = [get_devstatus(),
               "Environment :: Console",
               ## "Environment :: MacOS X",
               ## "Environment :: Win32 (MS Windows)",
               ## "Environment :: X11 Applications :: Qt",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Natural Language :: English",
               ## "Operating System :: Microsoft :: Windows",
               "Operating System :: POSIX :: Linux",
               ## "Operating System :: MacOS :: MacOS X",
               "Programming Language :: Python :: 2.7",
               ## "Programming Language :: Python :: 3.4",
               ## "Programming Language :: Python :: 3.5",
               "Topic :: Documentation :: Sphinx",
               "Topic :: Scientific/Engineering :: Physics",
               "Topic :: Software Development :: Libraries :: Python Modules"
               ]


##################
## Requirements ## 
##################
install_requires = ["numpy", "future", "scipy", "h5py", "fabio", "silx", "pyparsing", "PyMca5", "shapely", "matplotlib", "uncertainties","pint"]
extras_require = {"physics":["xraylib", "cctbx", "fdmnes"],\
                  "elastix":["SimpleITK"]\
                  }
setup_requires = ["testfixtures"]

###################
## Package setup ## 
###################
setup(name=PROJECT,
      version=get_version(),
      url="https://github.com/woutdenolf/spectrocrunch",
      author="Wout De Nolf",
      author_email="woutdenolf@users.sf.net",
      classifiers = classifiers,
      description="Spectroscopic imaging library (XRF/XAS)",
      long_description=get_readme(),
      install_requires=install_requires,
      extras_require=extras_require,
      setup_requires=setup_requires,
      packages=find_packages(),
      package_data={'spectrocrunch.resources': ['cif/*.*','id21/*.*']},
      license="MIT",
      cmdclass=cmdclass
      )
