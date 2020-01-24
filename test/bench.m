
bench_res := [* *];

for K in [GF(257), GF(521), GF(1031), GF(2053), GF(4099), GF(8209), GF(211^2), GF(521^2) ] do
    s := K.1;

    a := 2;
    R<x>:= PolynomialRing(K);

    h := x^5 + 1321*(s+112)*x^4 + 12321*(s+123)*x + 1765765 + s*1221;

    res := [* *];

    Append(~res, K);
    try
        time print (ZetaFunction(a, h)); // will include compile time
        print ("Compilation done");
        T := Cputime();
        ZetaFunction(a, h);
        t0 := Cputime(T);
        T := Cputime();
        ZetaFunction(a, h);
        t1 := Cputime(T);
        T := Cputime();
        ZetaFunction(a, h);
        t2 := Cputime(T);
        Append(~res, <h, a, (t0 + t1 + t2)/3, Max([t0, t1, t2]), Min([t0, t1, t2])>);
    catch e
        print (e);
    end try;

    h := x^7 + (s+7)*x^6 + 1321*(s+112)*x^5 + 12321*(s+123)*x + 1765765 + s*1221;
    try
        time print (ZetaFunction(a, h)); // will include compile time
        print ("Compilation done");
        T := Cputime();
        ZetaFunction(a, h);
        t0 := Cputime(T);
        T := Cputime();
        ZetaFunction(a, h);
        t1 := Cputime(T);
        T := Cputime();
        ZetaFunction(a, h);
        t2 := Cputime(T);
        Append(~res, <h, a, (t0 + t1 + t2)/3, Max([t0, t1, t2]), Min([t0, t1, t2])>);
    catch e
        print (e);
    end try;

    a := 3;
    try
        time print (ZetaFunction(a, h)); // will include compile time
        print ("Compilation done");
        T := Cputime();
        ZetaFunction(a, h);
        t0 := Cputime(T);
        T := Cputime();
        ZetaFunction(a, h);
        t1 := Cputime(T);
        T := Cputime();
        ZetaFunction(a, h);
        t2 := Cputime(T);
        Append(~res, <h, a, (t0 + t1 + t2)/3, Max([t0, t1, t2]), Min([t0, t1, t2])>);

    catch e
        print (e);
    end try;

    a := 4;
    try
        time print (ZetaFunction(a, h)); // will include compile time
        print ("Compilation done");
        T := Cputime();
        ZetaFunction(a, h);
        t0 := Cputime(T);
        T := Cputime();
        ZetaFunction(a, h);
        t1 := Cputime(T);
        T := Cputime();
        time ZetaFunction(a, h);
        t2 := Cputime(T);
        Append(~res, <h, a, (t0 + t1 + t2)/3, Max([t0, t1, t2]), Min([t0, t1, t2])>);

    catch e
        print (e);
    end try;
    
    Append(~bench_res, res);
end for;


/*
a:= 3;

R<x>:= PolynomialRing(GF(257));
h:= x^8 + x^7 + 1321*x^5 + 12321*x + 1765765;
time ZetaFunction(a, h);
print("ok");
time _ := ZetaFunction(a, h);
time _ := ZetaFunction(a, h);

R<x>:= PolynomialRing(GF(521));
h:= x^8 + x^7 + 1321*x^5 + 12321*x + 1765765;
time ZetaFunction(a, h);
time _ := ZetaFunction(a, h);
time _ := ZetaFunction(a, h);

R<x>:= PolynomialRing(GF(1031));
h:= x^8 + x^7 + 1321*x^5 + 12321*x + 1765765;
time ZetaFunction(a, h);
time _ := ZetaFunction(a, h);
time _ := ZetaFunction(a, h);

K:=GF(211^2);
b := K.1;
R<x>:= PolynomialRing(GF(211^2));
h:= x^7 + 1321*b*x^5 + (12321+b)*x + 1765765;
time ZetaFunction(a, h);
time _ := ZetaFunction(a, h);
time _ := ZetaFunction(a, h);

a:=2;
R<x>:= PolynomialRing(GF(2^32 + 15));
h := x^3 + 12321*x + 1765765;
time ZetaFunction(a, h);
print("ok");
time ZetaFunction(a, h);*/
