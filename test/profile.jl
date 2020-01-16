using Nemo, Coleman, Profile, ProfileView

a = 3
F, b = Nemo.FiniteField(211,2,"s")
R, x = PolynomialRing(F,"x")
h = x^7 + 1321*b*x^5 + (12321+b)*x + 1765765

Profile.clear();
@profile ZetaFunction(a, h);
println("Compilation done")
@time ZetaFunction(a, h)

Profile.clear();
@profile ZetaFunction(a, h);
ProfileView.view(C=true)
