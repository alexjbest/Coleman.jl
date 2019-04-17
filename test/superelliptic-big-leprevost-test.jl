using Nemo,Coleman


#p = 109
p = 35184372088891 # = 2^45 + 59

Br,w = FiniteField(p,1,"w")

Ri,x = PolynomialRing(PadicField(p, 3),'x')

#d = Int(degree(Br) - 1)
#bv = [w^i for i in 0:d]
#for xp in collect(Base.product([1:Int(characteristic(Br)) for i in 0:d]...))
#    a = sum([xp[i]*bv[i] for i in 1:d+1])
#    println(sqrtmod(f(a)), characteristic(Br))
#end

f=x^5 + 33//16*x^4 + 3//4*x^3  + 3//8*x^2 - 1//4*x + 1//16


P1=(-1,1)
P1n=(-1,-1)
P2=(0,ZZ(1)//4)
P2n=(0,ZZ(-1)//4)
P3=(1,-2)
P3n=(1,2)
pts = [P1,P1n,P2,P2n,P3, P3n]

@assert verify_pts(2, f, pts)

println(ColemanIntegrals(2, f, 2, p, 1, P1, P2))
