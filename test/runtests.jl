using Coleman
using Nemo
using Test
include("../src/LinearRecurrence.jl")

function test_zeta(a, h, inf_pts)
    print("zeta function y^",a,"=",h,"...")

    @test derivative(ZetaFunction(a, h))(ZZ(0)) == ZZ(count_points(a, h) + inf_pts)
    g = subst(h,gen(parent(h))+12)
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

    @test sort(rational_pts(a, h, bound)) == sort(rat_pts)

    println("PASS")
end

@testset "Coleman.jl" begin
    test_linear_recurrences_modx()

    R, x = PolynomialRing(FlintQQ, "x")
    # http://beta.lmfdb.org/Genus2Curve/Q/3153/a/9459/1
    h = x^5 +2*x^4 -3*x^3 -x^2 +x + 1//4
    test_rat_pts(2, h, 70, [(-1, -2 + 1//2), (-1, 1 + 1//2), (0, -1 + 1//2), (0, 0 + 1//2), (1, -1 + 1//2), (1, 0 + 1//2), (4, -37 + 1//2), (4, 36 + 1//2)])

    a = 2
    R, x = PolynomialRing(Nemo.GF(307), "x")

    h = x^9 + 1321*x^5 + 12321*x + 1765765

    test_zeta(a, h, 1)

    a = 3
    h = x^8 + 1321*x^5 + 12321*x + 1765765
    test_zeta(a, h, 1)

    R, x = PolynomialRing(Nemo.GF(601), "x")
    a = 4
    h = x^11 + 1765765
    test_zeta(a, h, 1)


end
