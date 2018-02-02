import spectrocrunch.math.linop as linop
import numpy as np
import sympy

ops = [linop.LinearOperator(1.3,0.1),linop.Clip(10,11),linop.Clip(3,7)]

x = -29
y = x
opc = linop.Identity()
for op in ops:
    y = op(y)
    opc = op*opc

print (opc)
print (y,opc(x))


