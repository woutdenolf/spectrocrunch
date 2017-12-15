import spectrocrunch.math.symbolic as symbolic

f = lambda x: 2*x
f = symbolic.Function(f,"f")

print(f)

f = 5*f

print(f)

f = f*f

print(f(5))
