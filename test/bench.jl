using Coleman
using Nemo

a = 3

R, x = PolynomialRing(Nemo.GF(2^32 + 15), "x")
h = x^5 + x^4 + 12321*x + 1765765
@time ZetaFunction(a, h) # will include compile time
println("Compilation done")
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)

R, x = PolynomialRing(Nemo.GF(257), "x")
h = x^8 + x^7 + 1321*x^5 + 12321*x + 1765765
@time ZetaFunction(a, h) # will include compile time
println("Compilation done")
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

#=
F, b = Nemo.FiniteField(211,2,"s")
R, x = PolynomialRing(F,"x")
h = x^8 + x^7 - b*1321*x^5 + b*12321*x + 1765765
@time ZetaFunction(a, h) # will include compile time
println("Compilation done")
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
=#

F = Nemo.GF(ZZ(2)^64 + 13)
R, x = PolynomialRing(F,"x")
h = x^4 + x^3 + 12321*x - 1765765
@time ZetaFunction(a, h) # will include compile time?
println("Compilation done")
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
