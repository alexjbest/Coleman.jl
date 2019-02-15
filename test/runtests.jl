using Coleman
using Nemo
using Test

function test_zeta(a, h, inf_pts)

    @test derivative(ZetaFunction(a, h))(ZZ(0)) == ZZ(count_points(a, h) + inf_pts)
    g = subst(h,gen(parent(h))+12)
    @test ZetaFunction(a, h) == ZetaFunction(a, g)

end

function test_rat_pts(a, h, bound, rat_pts)

    @test sort(rational_pts(a, h, bound)) == sort(rat_pts)

end

@testset "Coleman.jl" begin
    R, x = PolynomialRing(FlintQQ, "x")
    # http://beta.lmfdb.org/Genus2Curve/Q/3153/a/9459/1
    h =x^5 +2*x^4 -3*x^3 -x^2 +x + 1//4
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
