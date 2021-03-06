Change Log
==========

`develop`_ (unreleased)
-----------------------

* Noise propagation for fullfield XANES and micro-XANES (fluo and transmission)
* Alignment of XRF maps on raw spectrum (integrated with dead time correction, flux normalization and detector summing)
* HDF5 based job pipelines
* File system proxies


`0.0.1-alpha2`_ (unreleased)
----------------------------

* Cross-section calculations combining xraylib (tables) with fdmnes (calculation)
* Parse ESRF ID21 specific fullfield data
* fullfield XANES normalization
* Parse ESRF ID16b specific fluoXAS data
* ESRF ID21 specific data visualization
* 1D spectrum alignment (only phase correlation and simple properties)


`genesis`_ (unreleased)
-----------------------

* Initial version for Python 2.7
* Image alignment based on phase correlation, SIFT, Elastix and simple image properties (max, min, centroid and gaussfit)
* XRF fitting and deadtime correction of fluoXAS stacks using pyMca
* fluoXAS normalization
* Parse ESRF ID21 specific fluoXAS data


.. _genesis: https://github.com/woutdenolf/spectrocrunch/commit/genesis
.. _0.0.1-alpha2: https://github.com/woutdenolf/spectrocrunch/compare/genesis...v0.0.1-alpha2
.. _develop: https://github.com/woutdenolf/spectrocrunch/compare/v0.0.1-alpha2...master
