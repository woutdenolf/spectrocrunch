import spectrocrunch.materials.compoundfromname as compoundfromname
import spectrocrunch.materials.mixture as mixture
import spectrocrunch.materials.types as types

c1 = compoundfromname.compoundfromname("hematite")
c2 = compoundfromname.compoundfromname("calcite")
m = mixture.Mixture([c1,c2],[0.5,0.5],types.fractionType.weight)

m.markabsorber(symb="O",shells=["K","L2"],fluolines=["KA"])
m.markabsorber(symb="Fe",shells=["L1"])


print "\n".join(m.markinfo())

