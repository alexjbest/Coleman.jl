using Revise,Coleman,Nemo
set_printing_mode(FlintPadicField, :val_unit)
N = 5
a = 2
R,x = PolynomialRing(PadicField(67,N),"x")
Rb,xb = PolynomialRing(GF(67),"xb")
hb = xb^5 + 8
h = x^5 + 8
K = base_ring(h)
@assert N <= 8
#A = K(55 + 2*67 + 12*67^2 + 35*67^3 + 12*67^4 + 63*67^5 + 44*67^6 + 24*67^7)
#B = K(11 + 64*67 + 54*67^2 + 31*67^3 + 54*67^4 + 3*67^5 + 22*67^6 + 42*67^7 )
#E = (A,A)
#E2 = (B,B)
P = (K(1), K(3))

resd = rational_points(a, hb)

#Coleman.Chabauty(a,h,N - 1,Int(prime(K)),1,[([E , E2] , [:inf,:inf])],resd)
Coleman.Chabauty(a,h,N - 1,Int(prime(K)),1,[(P, :inf)],resd)
