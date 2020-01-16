


a:= 3;
/*
R<x>:= PolynomialRing(GF(257));
h:= x^8 + x^7 + 1321*x^5 + 12321*x + 1765765;
time ZetaFunction(a, h);
print("ok");
time ZetaFunction(a, h);
time ZetaFunction(a, h);
time ZetaFunction(a, h);

R<x>:= PolynomialRing(GF(521));
h:= x^8 + x^7 + 1321*x^5 + 12321*x + 1765765;
time ZetaFunction(a, h);
time ZetaFunction(a, h);
time ZetaFunction(a, h);

R<x>:= PolynomialRing(GF(1031));
h:= x^8 + x^7 + 1321*x^5 + 12321*x + 1765765;
time ZetaFunction(a, h);
time ZetaFunction(a, h);
time ZetaFunction(a, h);
*/

K:=GF(211^2);
b := K.1;
R<x>:= PolynomialRing(GF(211^2));
h:= x^7 + 1321*b*x^5 + (12321+b)*x + 1765765;
time ZetaFunction(a, h);
print("ok");
time ZetaFunction(a, h);
time ZetaFunction(a, h);
time ZetaFunction(a, h);

R<x>:= PolynomialRing(GF(2^32 + 15));
h := x^3 + 12321*x + 1765765;
time ZetaFunction(a, h);
print("ok");
time ZetaFunction(a, h);
