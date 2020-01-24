

for K in [GF(257), GF(521), GF(1031), GF(2053), GF(4099), GF(8209), GF(211^2), GF(521^2) ] do
    s := K.1;

    a := 2;
    R<x>:= PolynomialRing(K);

    h := x^5 + 1321*(s+112)*x^4 + 12321*(s+123)*x + 1765765 + s*1221;
    try
        time print (ZetaFunction(a, h)); // will include compile time
        print ("Compilation done");
        time ZetaFunction(a, h);
        time ZetaFunction(a, h);
        time ZetaFunction(a, h);
    catch e
        print (e);
    end try;

    h := x^7 + (s+7)*x^6 + 1321*(s+112)*x^5 + 12321*(s+123)*x + 1765765 + s*1221;
    try
        time print (ZetaFunction(a, h)); // will include compile time
        print ("Compilation done");
        time ZetaFunction(a, h);
        time ZetaFunction(a, h);
        time ZetaFunction(a, h);
    catch e
        print (e);
    end try;

    a := 3;
    try
        time print (ZetaFunction(a, h)); // will include compile time
        print ("Compilation done");
        time ZetaFunction(a, h);
        time ZetaFunction(a, h);
        time ZetaFunction(a, h);
    catch e
        print (e);
    end try;

    a := 4;
    try
        time print (ZetaFunction(a, h)); // will include compile time
        print ("Compilation done");
        time ZetaFunction(a, h);
        time ZetaFunction(a, h);
        time ZetaFunction(a, h);
    catch e
        print (e);
    end try;
end for;


/*
a:= 3;

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


K:=GF(211^2);
b := K.1;
R<x>:= PolynomialRing(GF(211^2));
h:= x^7 + 1321*b*x^5 + (12321+b)*x + 1765765;
time ZetaFunction(a, h);
print("ok");
time ZetaFunction(a, h);
time ZetaFunction(a, h);
time ZetaFunction(a, h);


a:=2;
R<x>:= PolynomialRing(GF(2^32 + 15));
h := x^3 + 12321*x + 1765765;
time ZetaFunction(a, h);
print("ok");
time ZetaFunction(a, h);*/
