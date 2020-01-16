using Coleman
using Hecke
using Nemo
using Test
include("../src/LinearRecurrence.jl")

function test_zeta(a, h, inf_pts)

    print("zeta function y^",a,"=",h,"...")

    Z = ZetaFunction(a, h)
    @info Z
    @info sqrt(ArbField(100)(order(base_ring(h))))

    @test IsWeil(numerator(Z), sqrt(ArbField(100)(order(base_ring(h)))))
    @test derivative(Z)(ZZ(0)) == ZZ(count_points(a, h) + inf_pts)
    g = subst(h,2*gen(parent(h))+12)
    @test ZetaFunction(a, h) == ZetaFunction(a, g)

    println("PASS")
end

function test_linear_recurrences_modx()
    print("linear recurrences mod x...")

    N = 3
    p = 7
    L = [11,13]
    R = [13,15]

    # Real base ring
    BrR= PadicField(p, N)#ResidueRing(ZZ,Int(p^N))
    RiR,xR = PolynomialRing(BrR,'x')

    FakeP,k = PolynomialRing(QQ, 'k')
    BrF= FmpqAbsSeriesRing(4,:k)
    k = gen(BrF)
    RiF,xF = PolynomialRing(BrF,'x')

    MF = matrix(RiF,2,2,[xF + k, xF, 1 + 3*k, xF + 1])
    MR = matrix(RiR,2,2,[xR + p, xR, 1 + 3*p, xR + 1])
    O = Coleman.LinearRecurrence(MF, L, R)

    o = [matrix(BrR,
                [ xadic_to_padic(t[i, j],BrR) for i = 1:nrows(t),j = 1:ncols(t)])
         for t in O]



    o == Coleman.LinearRecurrence(MR, L, R)
    println("PASS")
end

function test_rat_pts(a, h, bound, rat_pts)
    print("rational points y^",a,"=",h,"...")

    @test sort(rational_points(a, h, bound)) == sort(rat_pts)

    println("PASS")
end
function papertest()
    K = QadicField(41,2,20)
    R,x = PolynomialRing(K, "x")
    a = gen(K)
    h = x^4 + 7*x^3 + 3*x^2 - x
    P = Coleman.lift_point(3, h, K(1), 14*a+11) # the point (1,10^1/3)
    @test P[2]^3 == 10
    @test all([ColemanIntegrals(3, h, 3, 41, 2, P, :inf)[i] == 0 for i in RegularIndices(3, h)])
    println("PASS")
end

function rational_torsion_irrat_curve()
    A = 2
    K = QadicField(43,2,20)
    R,x = PolynomialRing(K, "x")
    D = Hecke.Hensel_factorization(x^2 -x+3)
    a = -coeff([D[k] for k in keys(D)][1],0)
    h = x^3 + (a-3)*x^2 + 16*a*x + 64
    P = (K(0),K(-8))
    @test all([ColemanIntegrals(A, h, 4, 43, 2, P, :inf)[i] == 0 for i in RegularIndices(A, h)])
end

function cubic_number_field_ex()
    A = 2
    K = QadicField(73,3,20)
    R,x = PolynomialRing(K, "x")
    D = Hecke.Hensel_factorization(x^3 -x^2+1)
    a = -coeff([D[k] for k in keys(D) if degree(D[k]) == 1][1],0)
    h = x^3 + (-3*a^2-2*a-3)*x^2 + (8*a^2+8*a-8)*x + (-16*a+16)
    P = (4*a^2, 4*a^2 - 4*a - 4)
    @test all([ColemanIntegrals(A, h, 4, 73, 3, P, :inf)[i] == 0 for i in RegularIndices(A, h)])
end

function quartic_number_field_ex()
    # LMFDB ecnf 4.4.725.1-89.1-a2
    A = 2
    K = QadicField(43,4,20)
    R,x = PolynomialRing(K, "x")
    D = Hecke.Hensel_factorization(x^4 - x^3 - 3*x^2 + x + 1)
    a = -coeff([D[k] for k in keys(D) if degree(D[k]) == 1][1],0)
    h =  x^3 + (2*a^3+6*a^2-9*a-4)*x^2 + (32*a^2-8)*x + (16*a^3+48*a^2-16*a-16)
    P = (K(0) , 4*a^2) # a 17-torsion point
    @test all([ColemanIntegrals(A, h, 4, 43, 4, P, :inf)[i] == 0 for i in RegularIndices(A, h)])
end


function bad()
    a = 3
    R, x = PolynomialRing(Nemo.FiniteField(67,2, "s")[1], "x")
    h = x^4 + 8
    test_zeta(a, h, 1)
end

function newbad()
    # Test 6
    a = 3
    K, s =Nemo.FiniteField(83,3,"s")
    R, x = PolynomialRing(K, "x")
    h = s*x^4 +x*3+ 8
    test_zeta(a, h, 1)
end

@testset "Coleman.jl" begin
    newbad()
    # Test 5
    bad()
    papertest()
    cubic_number_field_ex()
    rational_torsion_irrat_curve()

    test_linear_recurrences_modx()

    R, x = PolynomialRing(FlintQQ, "x")
    # http://beta.lmfdb.org/Genus2Curve/Q/3153/a/9459/1
    h = x^5 +2*x^4 -3*x^3 -x^2 +x + 1//4
    test_rat_pts(2, h, 70, [(-1, -2 + 1//2), (-1, 1 + 1//2), (0, -1 + 1//2), (0, 0 + 1//2), (1, -1 + 1//2), (1, 0 + 1//2), (4, -37 + 1//2), (4, 36 + 1//2)])

    # Test 1
    a = 3
    R, x = PolynomialRing(Nemo.GF(67), "x")
    h = x^5 + 8
    test_zeta(a, h, 1)

    # Test 2
    a = 2
    R, x = PolynomialRing(Nemo.GF(307), "x")
    h = x^9 + 1321*x^5 + 12321*x + 1765765
    test_zeta(a, h, 1)

    # Test 3
    a = 3
    h = x^8 + 1321*x^5 + 12321*x + 1765765
    test_zeta(a, h, 1)

    # Test 4
    R, x = PolynomialRing(Nemo.GF(601), "x")
    a = 4
    h = x^11 + 1765765
    test_zeta(a, h, 1)


end
