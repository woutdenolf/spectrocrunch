# -*- coding: utf-8 -*-

from . import compound


class CompoundFromList(compound.Compound):
    """Interface to a compound defined by a list of elements
    """

    def __init__(self, elements, frac, fractype, density=None, name=None):
        """
        Args:
            elements(list[str]): list of elements (["Fe","O"])
            frac(list[float]): element weight fractions or multiplicities
            ftype(types.fraction): element fraction type
            density(float): compound density (g/cm^3)
            name(Optional[str]): compound name
        """

        super(CompoundFromList, self).__init__(
            elements, frac, fractype, density=density, name=name
        )
