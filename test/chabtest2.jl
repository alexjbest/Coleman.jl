using Revise,Coleman,Nemo
set_printing_mode(FlintPadicField, :val_unit)
N = 3
a = 2
p = 31
R,x = PolynomialRing(PadicField(31,N),"x")
Rb,xb = PolynomialRing(GF(31),"xb")

hb = xb^5 + 12*xb^4 - 32*xb^3 - 512*xb^2 + 4096
h = x^5 + 12*x^4 - 32*x^3 - 512*x^2 + 4096
K = base_ring(h)
P = (K(0), K(64))

resd = rational_points(a, hb)

#Coleman.Chabauty(a,h,N - 1,Int(prime(K)),1,[([E , E2] , [:inf,:inf])],resd)
Coleman.Chabauty(a,h,N - 1,Int(prime(K)),1,[(P, :inf)],resd)
