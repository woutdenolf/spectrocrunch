# SpectroCrunch: spectroscopic imaging library (XRF/XAS)

> **Warning:**  
> This library is deprecated and no longer maintained.  
> It is replaced by [ewoksfluo](https://ewoksfluo.readthedocs.io) and [ewoksndreg](https://ewoksndreg.readthedocs.io).

## Getting started

### Install from PyPi

```bash
pip install [--user] spectrocrunch
```

### Install from source

```bash
git clone https://github.com/woutdenolf/spectrocrunch
cd spectrocrunch
pip install [--user] .
```

### Install non-pypi dependencies

```bash
git clone https://github.com/woutdenolf/spectrocrunch
. spectrocrunch/tools/linux-install-deps.sh [-u]
```

### Test

```bash
pip install -e .[test]
pytest .
```
