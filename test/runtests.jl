using Coleman
using Nemo
using Test

function test_zeta(a, h, inf_pts)

    @test derivative(ZetaFunction(a, h))(ZZ(0)) == ZZ(count_points(a, h) + inf_pts)

end

@testset "Coleman.jl" begin
    a = 2
    R, x = PolynomialRing(Nemo.GF(1009), "x")

    h = x^9 + 1321*x^5 + 12321*x + 1765765

    test_zeta(a, h, 1)

    a = 3
    h = x^8 + 1321*x^5 + 12321*x + 1765765
    test_zeta(a, h, 1)

    a = 4
    h = x^13 + 1765765
    test_zeta(a, h, 1)
end
