using Coleman
using Nemo

bench_res = []
for (K,s) in [(Nemo.GF(257),1), (Nemo.GF(521),1), (Nemo.GF(1031),1), (Nemo.GF(2053),1), (Nemo.GF(4099), 1), (Nemo.GF(8209), 1), Nemo.FiniteField(211,2,"s"), Nemo.FiniteField(521,2,"s") ]
    println("Doing $K")

    res = []

    push!(res, K)

    a = 2
    R, x = PolynomialRing(K, "x")

    h = x^5 + 1321*(s+112)*x^4 + 12321*(s+123)*x + 1765765 + s*1221
    try
        @time println(ZetaFunction(a, h)) # will include compile time
        println("Compilation done")
        t0 = @elapsed ZetaFunction(a, h)
        t1 = @elapsed ZetaFunction(a, h)
        t2 = @elapsed ZetaFunction(a, h)
        push!(res, (h, a, (t0 + t1 + t2)/3, max(t0, t1, t2), min(t0, t1, t2)))
    catch e
        println(e)
    end

    h = x^7 + (s+7)*x^6 + 1321*(s+112)*x^5 + 12321*(s+123)*x + 1765765 + s*1221
    try
        @time println(ZetaFunction(a, h)) # will include compile time
        println("Compilation done")
        t0 = @elapsed ZetaFunction(a, h)
        t1 = @elapsed ZetaFunction(a, h)
        t2 = @elapsed ZetaFunction(a, h)
        push!(res, (h, a, (t0 + t1 + t2)/3, max(t0, t1, t2), min(t0, t1, t2)))
    catch e
        println(e)
    end

    a = 3
    try
        @time println(ZetaFunction(a, h)) # will include compile time
        println("Compilation done")
        t0 = @elapsed ZetaFunction(a, h)
        t1 = @elapsed ZetaFunction(a, h)
        t2 = @elapsed ZetaFunction(a, h)
        push!(res, (h, a, (t0 + t1 + t2)/3, max(t0, t1, t2), min(t0, t1, t2)))
    catch e
        println(e)
    end

    a = 4
    try
        @time println(ZetaFunction(a, h)) # will include compile time
        println("Compilation done")
        t0 = @elapsed ZetaFunction(a, h)
        t1 = @elapsed ZetaFunction(a, h)
        t2 = @elapsed ZetaFunction(a, h)
        push!(res, (h, a, (t0 + t1 + t2)/3, max(t0, t1, t2), min(t0, t1, t2)))
    catch e
        println(e)
    end
    push!(bench_res, res)
end

bench_res

#=
a = 2
R, x = PolynomialRing(Nemo.GF(2^32 + 15), "x")
h = x^4 + 12321*x + 1765765
@time print(ZetaFunction(a, h)) # will include compile time
println("Compilation done")
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)

F = Nemo.GF(ZZ(2)^64 + 13)
R, x = PolynomialRing(F,"x")
h = x^4 + x^3 + 12321*x - 1765765
@time ZetaFunction(a, h) # will include compile time?
println("Compilation done")
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
@time ZetaFunction(a, h)
=#
