# SpectroCrunch: spectroscopic imaging library (XRF/XAS)

> **Warning:**  
> This library is deprecated and no longer maintained.  
> It is replaced by [ewoksfluo](https://ewoksfluo.readthedocs.io) and [ewoksndreg](https://ewoksndreg.readthedocs.io).

## Getting started

### Install dependencies:

```bash
git clone https://github.com/woutdenolf/spectrocrunch
. spectrocrunch/tools/linux-install-deps.sh [-u]
```

### Install from PyPi:

```bash
pip install [--user] spectrocrunch
```

### Install from source:

```bash
git clone https://github.com/woutdenolf/spectrocrunch
cd spectrocrunch
pip install [--user] .
```

### Test:

```bash
# Test source:
python setup.py test

# Test installation:
python -m spectrocrunch.tests.test_all

# Test with options:
python -m unittest -v spectrocrunch.tests.test_all.main_test_suite
pytest spectrocrunch
nose2 -v spectrocrunch

# Individual tests
python -m unittest spectrocrunch.align.tests.test_align.test_align.test_centroid
pytest spectrocrunch/align/tests/test_align.py::test_align::test_centroid
```

### Documentation:

- **Release:** [http://pythonhosted.org/spectrocrunch](http://pythonhosted.org/spectrocrunch)
- **Latest:** [http://spectrocrunch.readthedocs.io/en/latest/](http://spectrocrunch.readthedocs.io/en/latest/)

## Developers

[![Travis Status Master](https://travis-ci.org/woutdenolf/spectrocrunch.svg?branch=master)](https://travis-ci.org/woutdenolf/spectrocrunch?branch=master)  
[![Appveyor Status](https://ci.appveyor.com/api/projects/status/1txj75w5hjpmjfl3/branch/master?svg=true)](https://ci.appveyor.com/project/woutdenolf/spectrocrunch/branch/master)

- **Main development website:** [https://github.com/woutdenolf/spectrocrunch](https://github.com/woutdenolf/spectrocrunch)
- **Distribution website:** [https://pypi.python.org/pypi/spectrocrunch](https://pypi.python.org/pypi/spectrocrunch)
- **Developers guide:** [here](https://github.com/woutdenolf/wdncrunch/blob/master/tools/README.rst/)

## Use without installation

```bash
git clone https://github.com/woutdenolf/spectrocrunch
cd spectrocrunch
```

To import modules from a package without installing the package, add the directory of the package to the `PYTHONPATH` environment variable or add this to the top of your script:

```python
import sys
sys.path.insert(1, '.../spectrocrunch')
```

### Import as follows:

```python
from spectrocrunch.align.alignElastix import alignElastix as align
```

