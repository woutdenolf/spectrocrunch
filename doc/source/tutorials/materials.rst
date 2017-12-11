Calculating with materials
==========================

Materials can be defined by pure elements, chemical compounds and mixtures of compounds. 
    
.. code-block:: python

    import spectrocrunch.materials.element
    import spectrocrunch.materials.compoundfromformula
    import spectrocrunch.materials.mixture
    import spectrocrunch.materials.types
    import spectrocrunch.materials.plot

    element = spectrocrunch.materials.element.Element("Ca")
    
    compound1 = spectrocrunch.materials.compoundfromformula.CompoundFromFormula("PbSO4",density=6.29)
    
    compound2 = spectrocrunch.materials.compoundfromformula.CompoundFromFormula("CaSO4",density=2.32)
    
    mixture = spectrocrunch.materials.mixture.Mixture([compound1,compound2],[0.5,0.5],spectrocrunch.materials.types.fractionType.weight)


Cross sections
--------------

The main usage of these material classes is to calculate cross sections. The following cross-section are provided:
- Mass attenuation coefficient
- Mass absorption coefficient
- Partial photoionization cross section of a shell
- XRF cross section of an emission line
- Total scattering cross section
- Rayleigh cross section
- Compton cross section


XRF cross section
-----------------

We define the XRF cross section of one emission line of an element as the product of
- the partial photoionization cross section of the excited shell
- the fluorescence yield of the excited shell
- the radiative rate of the emission line
- the density of the material
- the mass fraction of the element


XRF spectra
-----------

.. code-block:: python

    import matplotlib.pyplot as plt

    fig,axs = plt.subplots(3,figsize=(6,9))
    
    plt.axes(axs[0])
    spectrum = element.lines(5,emin=2)
    spectrocrunch.materials.plot.xrf(spectrum)
    
    plt.axes(axs[1])
    spectrum = compound1.lines(5,emin=2)
    spectrocrunch.materials.plot.xrf(spectrum)
    
    plt.axes(axs[2])
    spectrum = mixture.lines(5,emin=2)
    spectrocrunch.materials.plot.xrf(spectrum)
    
    plt.tight_layout()
    plt.show()


