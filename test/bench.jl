using Coleman
using Nemo

a = 3


R, x = PolynomialRing(Nemo.GF(257), "x")
h = x^8 + x^7 + 1321*x^5 + 12321*x + 1765765
@time ZetaFunction(a, h) # will include compile time
print("Compilation done")
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)

R, x = PolynomialRing(Nemo.GF(521), "x")
h = x^8 + x^7 + 1321*x^5 + 12321*x + 1765765
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)

R, x = PolynomialRing(Nemo.GF(1031), "x")
h = x^8 + x^7 + 1321*x^5 + 12321*x + 1765765
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)

R, x = PolynomialRing(Nemo.FiniteField(211,2,"s")[1],"x")
h = x^8 + x^7 + 1321*x^5 + 12321*x + 1765765
@time ZetaFunction(a, h) # will include compile time
print("Compilation done")
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
