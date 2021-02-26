module Coleman

#  Copyright (C) 2009-2011 Moritz Minzlaff <minzlaff@daad-alumni.de>
#  Copyright (C) 2018-2019 Alex J. Best <alex.j.best@gmail.com>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#*****************************************************************************
#
# This module implements Minzlaff's algorithm for computing zeta functions of
# superelliptic curves in larger characteristic [1].
# This part of the code is due to Moritz Minzlaff.
# Porting of this code to Julia/Nemo, refactoring and modification to compute 
# Coleman primitives and Coleman integration functionality for superelliptic 
# curves is due to Alex Best.
#
# [1] Minzlaff, M.: Computing Zeta Functions of Superelliptic Curves in Larger
#   Characteristic. Mathematics in Computer Science. 3(2), 209--224 (2010)
#
#*****************************************************************************

include("LinearRecurrence.jl")
include("Misc.jl")

#import AbstractAlgebra.RelSeriesElem
#import AbstractAlgebra.Ring
#import AbstractAlgebra.PolyRing
#import AbstractAlgebra.Generic.LaurentSeriesFieldElem
using Hecke, Nemo, LoadFlint

const libflint = LoadFlint.libflint


export ColemanIntegrals, TinyColemanIntegralsOnBasis, ZetaFunction, AbsoluteFrobeniusActionOnLift, AbsoluteFrobeniusAction, TeichmullerPoint, lift_y, lift_x, verify_pt, verify_pts, count_points, rational_points, FrobeniusLift, padic_evaluate, BasisMonomials, RegularIndices, IsWeil

# Some dumb useless to everyone else functions that let me use nmod as if it were padic
function Nemo.frobenius(a::Union{Nemo.nmod, Generic.Res{fmpz}, SeriesElem, padic})
    return a
end
function Nemo.frobenius(a::Union{Nemo.nmod, Generic.Res{fmpz}, SeriesElem, padic}, b::Int)
    return a
end

function Nemo.frobenius(a::Generic.Poly)
    return parent(a)([frobenius(coeff(a, i)) for i in 0:degree(a)])
end

function Nemo.degree(R::ResRing{fmpz})
    return 1
end

function Nemo.degree(R::Nemo.NmodRing)
    return 1
end

#function Nemo.degree(R::Nemo.FlintPadicField)
#    return 1
#end

function (R::Nemo.FqFiniteField)(f::fmpz_poly)
    return f(gen(R))
end

# A few generalities on the differentials and the spaces W_{s,t}:
# The differential x^iy^j dx lies in
#    the iota-th block of W_{s,t}
# with
#    s \in \{ i, i-1, ..., max{-1, i- (b-1)} \}
#    t = (j div a)+1 if j \ge 0
#    t = (-j div a) if j < 0
#    iota = a - (j rem a) if j \ge 0
#    iota = (-j rem a) if j < 0

function Row(j,a)
#    Returns row index of x^...y^j dx
#
#    INPUT:
#
#    - ``j`` - a finite cardinal
#    - ``a`` - a finite cardinal

    if (j >= 0)
        return div(j, a)+1
    else
        return div((-j), a)
    end
end

function Block(j,a)
#    Returns block index of x^...y^j dx
#
#    INPUT:
#
#    - ``j`` - a finite cardinal
#    - ``a`` - a finite cardinal

    if (j >= 0)
        return a - mod(j, a)
    else
        return mod(-j, a)
    end
end

function ScalarCoefficients(j, k, a, hk, p, q, N)
# Returns the scalar coefficients \mu_{k,r,j}, i.e.
#    res_[r] = \mu_{k,r-1,j}
# (note the shift by 1 in the argument r!)
# the sequence has length b*k+1

    R1 = base_ring(hk)
    res_ = [zero(R1) for r in 1:degree(hk)+1]

    for r = 0:degree(hk)
        lambda = coeff(hk, r)
        # num = numerator of the two binomial expressions
        num = FlintQQ(1)
        for i = 0:N-2
            num = num*(-(j//a)-i)
        end
        # denom = denominator of the two binomial expressions
        denom = factorial(k)*factorial((N-1)-k)
        # last summand
        summand = num//denom
        sum = (-1)^(N-1+k)*summand
        # summing up going down
        for l = (N-2):-1:k
            summand = summand*(l+1-k)//(-(j//a)-l)
            sum = sum + (-1)^(l+k)*summand
        end
        sum = R1(numerator(sum))*inv(R1(denominator(sum)))
        res_[r+1] = p*lambda*sum
    end

    return res_
end

function myinv(M)
    de = det_df(M)
    ad = zero(M)

    for i in 1:nrows(M)
        for j in 1:ncols(M)
            ad[i,j] = (-1)^(i+j)*det_df(M[[q for q in 1:nrows(M) if q != i], [w for w in 1:nrows(M) if w != j]])
        end
    end

    @assert inv(de)*transpose(ad)*M == 1
    return inv(de)*transpose(ad)
end

function RSCombination(h)
# returns polynomial sequences r_ and s_ such that
# r_[i]h + s_[i]Derivative(h) = x^i
# where i \le b-2 and deg r_[i] \le b-2. deg s_[i] \le b-1

    b = degree(h)
    rk = b+(b-1)
    R = base_ring(parent(h))
    RMat = MatrixSpace(R, rk, rk)

    dh = derivative(h)

    M = zero(RMat)
    for c = 1:rk
        for r = 1:b-1
            if ((c-1)-(r-1) >= 0)
                M[r,c] = coeff(h, (c-1)-(r-1))
            end
        end
        for r = b:rk
            if ((c-1)-(r-b) >= 0)
                M[r,c] = coeff(dh, (c-1)-(r-b))
            end
        end
    end
    try
        Mi = myinv(M)
        #@info ">>>>>>>>>>>>>>2"
        #@info M
        #@info Mi
    catch e
        _Mi, _d = pseudo_inv(lift(M))
        Mi = inv(R(_d)) * RMat(_Mi)
    end
    if Mi isa Tuple
        Mi = divexact(Mi[1],Mi[2]) # TODO be careful here with p-adic prec?
    end

    resR_ = [ parent(h)([ Mi[i,c] for c in 1:(b-1) ]) for i in 1:(b-1) ]
    resS_ = [ parent(h)([ Mi[i,c] for c in b:rk ]) for i in 1:(b-1) ]

    return resR_, resS_
end

function CastElement(R, e)
#    Return the image / a preimage of e in R
#
#    INPUT:
#
#    -  ``R`` - (a polynomial ring over) UnramifiedQuotientRing(K,N)
#    -  ``e`` - an element of (or a polynomial over)
#               UnramifiedQuotientRing(K,N') with N' = N-1 or N+1
#
#    OUTPUT:
#
#    An element of R
#
#    NOTE:
#
#    This function is needed since Magma cannot coerce between
#    UnramifiedQuotientRing(K,N) with different N

    RR = base_ring(R)
    # If K/L with L prime, then RR is the UnramifiedQuotientRing(L,N)
    # If R is not a polynomial ring, then (RR eq RRR)
    RRR = base_ring(RR)
    e_ = Eltseq(e)
    res_ = []
    for i = 1:length(e_)
        e__ = Eltseq(e_[i])
        res__ = []
        for j = 1:length(e__)
            res__[j] = RRR(e__[j])
        end
        res_[i] = RR(res__)
    end
    return R(res_)
end


function CastMatrix(R, M)
#    Return the image / a preimage of M over R
#
#    INPUT:
#
#    -  ``R`` - a (polynomial) matrix ring over UnramifiedQuotientRing(K,N)
#    -  ``M`` - a (polynomial) matrix over UnramifiedQuotientRing(K,N') with
#               N' = N-1 or N+1
#
#    OUTPUT:
#
#    An element of R
#
#    NOTE:
#
#    This function is needed since Magma cannot coerce between
#    UnramifiedQuotientRing(K,N) with different N
#
    RR = base_ring(R)
    res = zero_matrix(RR, nrows(M), ncols(M))
    for i = 1:nrows(M)
        for j = 1:ncols(M)
            res[i,j] = cast_poly_nmod(RR, M[i,j])
        end
    end
    return res
end

function CastBaseMatrix(R, M)
#    Return the image / a preimage of M over R
#
#    INPUT:
#
#    -  ``R`` - a (polynomial) matrix ring over UnramifiedQuotientRing(K,N)
#    -  ``M`` - a (polynomial) matrix over UnramifiedQuotientRing(K,N') with
#               N' = N-1 or N+1
#
#    OUTPUT:
#
#    An element of R
#
#    NOTE:
#
#    This function is needed since Magma cannot coerce between
#    UnramifiedQuotientRing(K,N) with different N
#
    RR = base_ring(R)
    res = zero_matrix(RR, nrows(M), ncols(M))
    for i = 1:nrows(M)
        for j = 1:ncols(M)
            res[i,j] = RR(lift_elem(M[i,j]))
        end
    end
    return res
end

function HRedMatrix(t, iota, a, h, R1PolMatH, pts)
# given row index t and block index iota,
# the equation of the curve (via a,h)
# return the horizontal reduction matrix
# for row t and block iota
# also return the denominators as a sequence of polynomials
# i.e.
# resM = M_H^{t,\iota}(s)
# resD = d_H^{t,\iota}(s)
#
    R1Pol = parent(h)
    s = gen(R1Pol)

    resM = zero(R1PolMatH)

    b = degree(h)
    lambda = lead(h)
    h1 = h - lambda*s^b # h - leading term (h)

    resD =  lambda*(b*(a*t+iota-a) -a*s)
    c_ = [ a*coeff(h1, 0)*s ]
    c_ = vcat(c_, [ R1Pol(a*coeff(h1, i)*s -
                          (a*t+iota-a)*coeff(derivative(h1), i-1)) for i in 1:(b-1) ])

    for i = 1:b-1
        resM[i,i+1] = resD
    end
    for i = 1:b
        resM[b,i] = c_[i]
    end
    for i = 1:length(pts)
        resM[b,b+i] = -a*(pts[i][2])^(-(iota-a))
        resM[b+i, b+i] = (resD)*pts[i][1]
    end

    return resM, resD
end

function HRedMatrixSeq(genM, genD, L_, R_, DDi, slr, p, N, B, Vi, R1MatH,
R0PolMatH)
#    Given the generic reduction matrix genM for the current row and block
#    together with its denominator genD
#    and given interval boundaries L_ and R_
#    return the matrix sequences specified by these intervals, i.e.
#    resM_[l] = M_H^{t,\iota}(l) and resD_[l] = "d_H^{t,\iota}(l)"
#
#    NOTE:
#
#    Computations are carried out mod p^N
#    but the result is given mod p^{N+1}
#
    R1Pol = parent(genD)
    R1 = base_ring(R1Pol)

    R0Pol = base_ring(R0PolMatH)
    R0 = base_ring(R0Pol)
    R0PolMat = MatrixSpace(R0Pol, 1, 1)

    tempM_ = LinearRecurrence(transpose(CastMatrix(R0PolMatH,genM)), L_,
                              R_, DDi, slr)
    tempM_ = [ transpose(tempM_[m]) for m in 1:length(tempM_) ]
    tempD_ = LinearRecurrence(transpose(R0PolMat(cast_poly_nmod(R0Pol,genD))),
                               L_, R_, DDi, slr)
    tempD_ = [ transpose(tempD_[m]) for m in 1:length(tempD_) ]
    if (N < B)    # we need to compute the remaining matrices
        if (N == 1)    # everything is congruent mod p
            tempM_ = vcat(tempM_, [ tempM_[1] for l in (N+1):B ])
            tempD_ = vcat(tempD_, [ tempD_[1] for l in (N+1):B ])
        else    # apply the vandermonde trick
                # denominators
            R0Mat = parent(tempD_[1])
            tempD_ = vcat(tempD_, [ zero(R0Mat) for l in (N+1):B ])
            taylor_ = [zero(R0Mat) for l in 1:N]
            for l = 1:N
                for m = 1:N
                    taylor_[l] = taylor_[l] + tempD_[m]*Vi[m,l]
                end
            end
            for l = N+1:B
                tempD_[l] = zero(R0Mat)
                c = one(R0)
                for i = 1:N
                    tempD_[l] = tempD_[l] + taylor_[i]*c
                    c = c*FlintZZ(l) # Ideally we should be able to write c*l here?
                end
            end
            # matrix
            R0Mat = parent(tempM_[1])
            tempM_ = vcat(tempM_, [ zero(R0Mat) for l in (N+1):B ])
            taylor_ = [zero(R0Mat) for l in 1:N]
            for l = 1:N
                for m = 1:N
                    taylor_[l] = taylor_[l] + tempM_[m]*Vi[m,l]
                end
            end
            for l = N+1:B
                tempM_[l] = zero(R0Mat)
                c = one(R0)
                for i = 1:N
                    tempM_[l] = tempM_[l] + taylor_[i]*c
                    c = c*FlintZZ(l) # Ideally we should be able to write c*l here?
                end
            end
        end
    end
    resM_ = [ CastBaseMatrix(R1MatH,tempM_[l]) for l in 1:B ]
    resD_ = [ R1(lift_elem(tempD_[l][1,1])) for l in 1:B ]

    return resM_, resD_
end

function HReduce(i, b, iota, mu_, genM, genD, M_, D_, p, R1ModH)
    # reduces the differential T_{(i,j),k} horizontally
    #
    R1 = base_ring(R1ModH)

    res = zero(R1ModH)

    # Note: #mu_ = b*k+1
    res[1,1] = mu_[end]   # Recall: mu_[r] = mu_{k,r+1,j}
    ##@info "res"
    ##@info res

    for l = (i+ length(mu_)):-1:1
        for m = 1:b-1
            res = mul!(res, res, Evaluate(genM, R1(p*l-m)))
            res *= inv(evaluate(genD, R1(p*l-m)))
        end
        ##@info "res",res;
        res *= Evaluate(genM, R1(p*l-b))
        ##@info "res",res;
        d = evaluate(genD, R1(p*l-b))
        ##@info "d",d;
        res = R1ModH([ R1(lift_elem(divexact(res[1,m],d))) for m in 1:R1ModH.ncols ])
        ##@info "res",res;
        res *= M_[l]
        res *= inv(D_[l])
        res = mul!(res, res, Evaluate(genM, R1((l-1)*p)))
        res *= inv(evaluate(genD,R1((l-1)*p)))
        if ((l-1)-i-1 >= 0)
            res[1,1] += mu_[(l-1)-i]
        end
    end

    return res
end

function VRedMatrixSeq(j, a, h, r_, s_, p, N, R1MatV, R1PolMatV, pts)
# Given the data to compute the generic reduction matrix
# (and its denominator) of the iota-th block,
# return the matrix sequences needed for vertical reduction, i.e.
## resM_[k] = M_V^{\iota}(k) and resD_[k] = "d_V^{\iota}(k)"
#
    b = degree(h)
    R1 = base_ring(h)
    R1Pol = parent(h)
    t = gen(R1Pol)
    R1PolMat = MatrixSpace(R1Pol, 1, 1)

    t_ = vcat([ 0 ], [ Row(-p*(a*k + j), a) for k in 0:(N-1) ])
    L_ = [ t_[i] for i in 1:(length(t_)-1) ]
    R_ = [ t_[i] for i in 2:length(t_) ]
    slr = floor(Int, log(4, R_[end]))
    DDi = UpperCaseDD(one(R1), R1(2^slr), 2^slr)
    DDi = inv(DDi)
    #@info "DDi"
    #@info DDi

    iota = Block(-p*j, a)

    genM = zero(R1PolMatV)
    for i = 1:b-1
        for m = 1:b-1
            genM[i,m] = (a*t + iota-a)*coeff(r_[i], m-1) +
                  a*coeff(derivative(s_[i]), m-1)
              #@info genM, coeff(r_[i], m-1), coeff(derivative(s_[i]), m-1)
        end
    end
    for m = 1:length(pts)
        for i = 1:b-1
            genM[i,m+b-1] = -a*evaluate(s_[i], pts[m][1])*(pts[m][2])^(-iota)
        end

        genM[b-1+m, b-1+m] = (a*t +iota-a)*(pts[m][2])^(-a)

    end
    #@info "genM"
    #@info genM
    #@info L_
    #@info R_
    #@info DDi
    #@info slr
    resM_ = LinearRecurrence(transpose(genM), L_, R_, DDi, slr)
    #@info resM_
    map!(transpose, resM_,resM_)

    genD = R1PolMat(a*t +iota-a)
    #@info "genD"
    #@info R_
    #@info L_
    #@info genD
    tempD_ = LinearRecurrence(transpose(genD), L_, R_, DDi, slr)
    #@info tempD_
    map!(transpose, tempD_, tempD_)
    #@info tempD_
    resD_ = [ tempD_[k][1,1] for k in 1:N ]
    #@info resM_
    #@info resD_

    return resM_, resD_
end

_lift_to(::Type{padic}) = fmpq

_lift_to(::Type{qadic}) = fmpz_poly

function VReduce(i, j, a, h, wH_, M_, D_, R1ModV)
    # "vertically" reduces the already
    # "horziontally reduced" differential
    # w_{(i,j)} = wH_[*,j,i+1]
    #
    R1 = base_ring(h)

    b = degree(h)
    N = length(wH_)

    res = R1ModV([ wH_[N][j][i+1][1,m] for m in 2:R1ModV.ncols+1 ])
    #@info res

    tmp = _lift_to(elem_type(R1))()

    for k in (N-1):-1:1
        res = mul!(res, res, M_[k+1])
        d = D_[k+1]
        ##@info d
        
        t = Vector{elem_type(R1)}(undef, R1ModV.ncols)
        #@show @which lift_elem(divexact(res[1, m], d))
        for m in 1:R1ModV.ncols
          z = lift_elem!(tmp, divexact(res[1, m], d))
          t[m] = R1(z)
        end
        #res = R1ModV([ R1(lift_elem(divexact(res[1,m], d))) for m in 1:R1ModV.ncols ])
        res = R1ModV(t)
        ##@info res
        # Add new term
        res = R1ModV([ wH_[k][j][i+1][1,m] + res[1,m-1] for m in 2:R1ModV.ncols+1 ])
        ##@info res
    end

    res *= M_[1]
    #@info "M_",M_[1]
    res *= inv(D_[1])
    #@info "D_",D_[1]

    return res
end

function lift_fq_to_qadic(R, a)
    if typeof(a) <: Union{<: ResElem, Nemo.gfp_elem}
        return R(lift_elem(a))
    else
        t = FmpzPolyRing(FlintZZ,:x)([coeff(a, i) for i in 0:degree(R)-1])
        if degree(R) == 1
            return R(coeff(t,0))
        end
        return R(t)
    end
end

function lift_fq_to_qadic_poly(R::PolyRing, f)
    #Ry, _ = PolynomialRing(ResidueRing(FlintZZ, characteristic(base_ring(parent(f)))^N), "y")
    return R([lift_fq_to_qadic(base_ring(R), coeff(f, i)) for i in 0:degree(f)])
end

function AbsoluteFrobeniusAction(a, hbar, N, pts = [])
    K = base_ring(hbar)
    p = convert(Int,characteristic(K))
    n = degree(K)

    if n == 1
        if fits(Int, FlintZZ(p)^(N+1))
            R0 = ResidueRing(FlintZZ, p^N)
            R1 = ResidueRing(FlintZZ, p^(N+1))
        else
            R0 = ResidueRing(FlintZZ, FlintZZ(p)^N)
            R1 = ResidueRing(FlintZZ, FlintZZ(p)^(N+1))
        end
    else
        R0 = FlintQadicField(p, n, N)
        R1 = FlintQadicField(p, n, N + 1)
    end
    R0Pol,t1 = PolynomialRing(R0,'t')
    R1Pol,t2 = PolynomialRing(R1,'t')


    h = lift_fq_to_qadic_poly(R1Pol, hbar)
    #@info h

    return AbsoluteFrobeniusActionOnLift(a, h, N, p, n, pts)
end

function AbsoluteFrobeniusActionOnLift(a, h, N, p, n, pts = [])#(a::RngIntElt, hbar::RngUPolElt,N::RngIntElt)\
#-> AlgMatElt
#
#   Implements [1, Algorithm 1]

#   INPUT:

#   -  ``a`` - an integer > 1
#   -  ``hbar`` - a squarefree univariate polynomial over a finite field
#                 of degree coprime to a
#   -  ``N`` - an integer > 0 setting the desired precision

#   OUTPUT:

#   A integer matrix modulo p^N representing the action of the
#   absolute Frobenius on the first crystalline cohomology space
#   of the smooth projective model of y^a - hbar = 0.

#   NOTE:

#   The complexity is O( p^(1/2) n MM(g) N^(5/2) + \log(p)ng^4N^4 )
#
    # Step 0: Setup
    b = degree(h)
    l = length(pts)

    # Check user input
    #(! IsFinite(K)) && error("The curve must be defined over a finite field.")
    #(! IsSeparable(h)) && error("The current implementation only supports squarefree h.")
    (gcd(a,b) != 1) && error("The current implementation needs $a and the degree of $h to be coprime.")
    (a < 2) && error("Please enter an integer a > 1.")
    (b < 2) && error("Please enter a polynomial h of degree > 1.")
    (N < 1) && error("Please enter a positive precision N")
    if p isa Union{Integer,fmpz}
        q = p^n
        (p <= (a*N-1)*b) && error("Characteristic too small", (a*N - 1)*b)

        if n == 1
            if true
                R0 = FlintPadicField(p, N)
                R1 = FlintPadicField(p, N + 1)
            elseif fits(Int, FlintZZ(p)^(N+1)) # Old code, maybe more efficient eventually
                R0 = ResidueRing(FlintZZ, p^N)
                R1 = ResidueRing(FlintZZ, p^(N+1))
            else
                R0 = ResidueRing(FlintZZ, FlintZZ(p)^N)
                R1 = ResidueRing(FlintZZ, FlintZZ(p)^(N+1))
            end
        else
            R0 = FlintQadicField(p, n, N)
            R1 = FlintQadicField(p, n, N + 1)
        end
    else # use a power series ring!
        if n == 1
            R0,_ = PowerSeriesRing(FlintQQ, N,     "p", cached=true, model=:capped_absolute)
            R1,_ = PowerSeriesRing(FlintQQ, N + 1, "p", cached=true, model=:capped_absolute)
            q = p
        else
            error("No q-adic p as a variable yet")
        end
    end

    R0Pol,t1 = PolynomialRing(R0,'t')
    R1Pol,t2 = PolynomialRing(R1,'t')

    Rt,t3 = PolynomialRing(FlintZZ,'t')
    h = cast_poly_nmod(R1Pol,h)
    #@info h

    pts = [(R1(lift_elem(P[1])),R1(lift_elem(P[2]))) for P in pts]

    # Step 1: Horizontal reduction
    R1MatH = MatrixSpace(R1, b + l, b + l)
    R1ModH = MatrixSpace(R1, 1, b + l)
    R1PolMatH = MatrixSpace(R1Pol, b + l, b + l)
    R0PolMatH = MatrixSpace(R0Pol, b + l, b + l)

    wH_ = [ [ [] for j in 1:(a-1) ] for k in 0:(N-1) ]
    # stores the results of the reduction
    # wH_[k+1,j,i+1] = w_{(i,j),k}
    # Note: w_{(i,j),k} is nonzero only in the
    # iota(j)-th block, so _only_ this block is stored

    # vandermonde trick: preliminaries
    R0Mat = MatrixSpace(R0, N, N)
    if (N < b-1 +b*(N-1)) && (N > 1)
        V = R0Mat( [ i^j for j in 0:N-1 for i in 1:N ])
        Vi = inv(V)
        if Vi isa Tuple
            Vi = divexact(Vi[1],Vi[2])
        end
    else
        Vi = one(R0Mat)
    end

    hk = one(R1Pol)
    hFrob = R1Pol([ frobenius(coeff(h,i)) for i in 0:degree(h) ])
    # at the start of the k-th loop hk = (hFrob)^k
    for k = 0:(N-1)
        # reduction matrix sequences: preliminaries
        B = b-1 +b*k
        mn = min(N, B)
        L_ = [ (l-1)*p for l in 1:mn ]
        R_ = [ (l*p -b-1) for l in 1:mn ]
        slr = floor(Int, log(4, R_[end]))
        DDi = UpperCaseDD(one(R0), R0(2^slr), 2^slr)
        DDi = inv(DDi)
        #@info "DDi"
        #@info DDi

        for j = 1:a-1
            # j and k fix the row index
            t = Row(-p*(a*k +j), a)
            # horizontal reductions are performed
            # row by row from "bottom to top"

            #iota = Block(-p*(a*k +j), a)
            # j fixes the block index
            # Note: this really is independent of k!
            iota = Block(-p*j, a)
            # Block(-p*(a*k+j), a) = Block(-p*j, a)
            @assert( -(t*a+iota) == -p*(a*k+j) )

            # generic reduction matrix
            genM, genD = HRedMatrix(t, iota, a, h, R1PolMatH, pts)
            #@info "genM"
            #@info genM,genD

            # reduction matrix sequences: evaluation
            M_, D_ = HRedMatrixSeq(genM, genD, L_, R_, DDi, slr,
                                   p, N, B, Vi, R1MatH, R0PolMatH)
            ##@info "M_"
            ##@info M_,D_

            # approximate frobenius action
            mu_ = ScalarCoefficients(j, k, a, hk, p, q, N)
            #@info "Mu_"
            #@info mu_

            # reduce
            wH_[k+1][j] = [ HReduce(i, b, iota, mu_, genM, genD, M_,
                                    D_, p, R1ModH) for i in 0:(b-2) ]
        end
        hk *= hFrob
    end
    #@info "wH_"
    #@info wH_


    # Step 2: Vertical reduction
    R1MatV = MatrixSpace(R1, b-1 + l, b-1 + l)
    R1ModV = MatrixSpace(R1, 1, b-1 + l)
    R1PolMatV = MatrixSpace(R1Pol, b-1 + l, b-1 + l)
    wV_ = [ [] for j in 1:(a-1) ]
    # stores the results of the reduction
    # wV_[j,i+1] = w_{(i,j)}
    # Note: w_{(i,j)} is nonzero only in the
    # iota(j)-th block, so _only_ this block is stored
    # Note: block size is now b-1!
    # (as opposed to b during horizontal reduction)

    # reduction matrix sequences: preliminaries
    # compute the r_i and s_i needed to define the
    # vertical reduction matrices
    r_, s_ = RSCombination(h)
    #@info "RS"
    #@info r_,s_

    for j = 1:a-1
        # reduction matrix sequences: evaluation
        M_, D_ = VRedMatrixSeq(j, a, h, r_, s_, p, N,
                                R1MatV, R1PolMatV, pts)
        #@info "MD"
        #@info M_
        #@info "MD"
        #@info D_

        # reduce
        wV_[j] = [ VReduce(i, j, a, h, wH_, M_, D_, R1ModV) for i in 0:(b-2) ]
    end
    #@info "wV_", wV_

    # Step 3: Assemble output
    R0Mat = MatrixSpace(R0, (a-1)*(b-1), (a-1)*(b-1))
    res = zero(R0Mat)
    for j = 1:a-1
        for i = 0:b-2
            for m = 1:b-1
                res[((j-1)*(b-1) +i+1), ((Block(-p*j, a)-1)*(b-1) +m)] = R0(lift_elem(wV_[j][i+1][1,m]))
            end
        end
    end
    #@info "res"
    #@info res
    @assert res != 0

    # Return just the matrix of frobenius if we have no points
    if l == 0
        return res
    end

    # Get the evaluations

    R0ColMat = MatrixSpace(R0, (a-1)*(b-1), l)
    col = zero(R0ColMat)
    for j = 1:a-1
        for i = 0:b-2
            for m = 1:l
                col[((j-1)*(b-1) +i+1), m] = R0(lift_elem(wV_[j][i+1][1,(b-1+m)]))
            end
        end
    end

    return res,col
end

function ZetaFunction(a, hbar)#(a::RngIntElt, hbar::RngUPolElt)
#
#   Implements [1, Corollary]

#   INPUT:

#   -  ``a`` - an integer > 1
#   -  ``hbar`` - a squarefree univariate polynomial over a
#                 finite field of degree coprime to a

#   OUTPUT:

#   A rational function over FlintQQ
#
    (gcd(a,degree(hbar)) != 1) && error("The current implementation needs $a and the degree of $hbar to be coprime.")
    # Step 0: Setup
    p = convert(Int,characteristic(base_ring(hbar)))
    q = order(base_ring(hbar))
    n = degree(base_ring(hbar))
    g = divexact(((a-1)*(degree(hbar)-1)),2)

    # Step 1: Determine needed precision
    bound = n*g/2 + 2*g*log(p,2)
    # N is the first integer strictly larger than bound
    N = floor(Int, bound+1)
    #@info N

    # Step 2: Determine absolute Frobenius action mod precision
    M = AbsoluteFrobeniusAction(a, hbar, N)

    #@info M
    # Step 3: Determine Frobenius action mod precision
    MM = deepcopy(M)
    for i in 1:n-1
        # Apply Frobenius to MM
        map!(frobenius, MM, MM)
        # Multiply
        M = MM * M
    end
    #@info M

    # Step 4: Determine L polynomial
    ZPol,t = PolynomialRing(FlintZZ,"t")
    #CP = charpoly(PolynomialRing(base_ring(M),"t")[1],M::MatElem{RingElem})
    #CP = invoke(charpoly, Tuple{Ring, Union{MatElem{Nemo.nmod},Generic.Mat}},  PolynomialRing(base_ring(M),"t")[1], M)
    CP = charpoly(PolynomialRing(base_ring(M),"t")[1], M)
    #@info CP
    Chi = cast_poly_nmod(ZPol, CP)
    L = numerator(t^(2*g)*(Chi)(1//t))
    coeff_ = [ coeff(L, i) for i in 0:(2*g) ]
    prec = FlintZZ(p)^N
    mid = prec >> 1
    for i = 0:g
        if (coeff_[i+1] > mid)
            coeff_[i+1] = coeff_[i+1]-prec
        end
    end
    for i = 0:g-1
        coeff_[2*g-i+1] = (q^(g-i))*coeff_[i+1]
    end
    L = ZPol(coeff_)

    # Step 5: Output zeta function
    return L // (q*t^2 - (q+1)*t + 1)
end



function IsWeil(P, sqrtq)
    (discriminant(P) == 0) && error("Polynomial not squarefree, so root-finding is hard?")


    prec = 100
    rts = []
    while true
        R = AcbPolyRing(AcbField(prec),:x)
        try
            rts = roots(R(P))
            break
        catch e
            prec *= 2
        end
    end
    Q = AcbField(prec)

    return all([overlaps(abs(a),abs(Q(sqrtq)^(-1))) for a in rts])

end

function Nemo.root(a::Union{Nemo.padic,Nemo.qadic}, n::Int)
    try
        if valuation(a) != 0
            if valuation(a) % n == 0
                K = parent(a)
                N = divexact(valuation(a), n)
                uniformizer = K(prime(K))
                return exp(log(a//(uniformizer^valuation(a)))//n)*uniformizer^N
            else
                error("no root over ground field")
            end
        end
        return exp(log(a)//n)
    catch
        K = parent(a)
        R,x = PolynomialRing(K, "x")
        D =  Hecke.Hensel_factorization(x^n-a)
        return -coeff([D[k] for k in keys(D) if degree(D[k]) == 1][1],0)
    end
end

function Nemo.root(r::Nemo.gfp_elem, a::Int)
    K = parent(r)
    for x in 0:Int(characteristic(K))-1
        if K(x)^a == r
            return K(x)
        end
    end
    error("no root")
end

function Nemo.root(r::Nemo.FinFieldElem, a::Int; not = nothing)
    K = parent(r)
    u = gen(K)
    p = characteristic(K)
    n = degree(K)
    for x in K
        if x^a == r && (not == nothing || x != not)
            return x
        end
    end
    error("no root")
end

function count_points_bad(a::Int, h::PolyElem{Nemo.gfp_elem})
    return length(rational_points(a, h))
end

function count_points(a::Int, h::PolyElem{Nemo.gfp_elem})
    K = base_ring(h)
    N = gcd(Int(characteristic(K)) - 1, a)
    su = 0
    for x in 0:Int(characteristic(K))-1
        try
            y = root(h(x), a)
            if y == 0
                su += 1
            else
                su += N
            end
        catch
        end
    end
    return su
end

function count_points(a::Int, h::PolyElem{<:Nemo.FinFieldElem})
    K = base_ring(h)
    N = gcd(Int(order(K)) - 1, a)
    su = 0
    for x in K
        try
            y = root(h(x), a)
            if y == 0
                su += 1
            else
                su += N
            end
        catch
        end
    end
    return su
end

function rational_points(a::Int, h::PolyElem{Nemo.gfp_elem})
    K = base_ring(h)
    N = gcd(Int(characteristic(K)) - 1, a)
    zeta = gen(K)^Int(divexact(characteristic(K) - 1, N))
    ret = []
    for x in 0:Int(characteristic(K))-1
        try
            y = root(h(x), a)
            if y == 0
                push!(ret, (x, y))
            else
                push!(ret, [(x,zeta^i*y) for i in 0:N-1]...)
            end
        catch
        end
    end
    return ret
end


function superelliptic_automorphism(a, h, p, n, P)
    @assert mod(ZZ(p^n - 1),ZZ(a)) == 0
    K = base_ring(h)
    R,X = PolynomialRing(K, "x")
    D = Hecke.Hensel_factorization(divexact(X^a - 1, X - 1))
    #@info D
    zeta = -coeff([D[k] for k in keys(D) if degree(D[k]) == 1][1],0)

    return (P[1], zeta*P[2])
end

function FrobeniusLift(a::Int, h, p::Int, P::Tuple{Any,Any}, m = 1)
    if m != 1 # apply recursively
        return FrobeniusLift(a, h, p, FrobeniusLift(a, h, p, P, m - 1), 1)
    end
    K = base_ring(h)
    P = (K(P[1]), K(P[2]))
    return (frobenius((P[1])^p, degree(K) - 1),
            frobenius(P[2]^p * root(1 + (frobenius(h)(P[1]^p) - h(P[1])^p)//(h(P[1])^p), a), degree(K) - 1))

end

function TeichmullerPoint(a::Int, h, p::Int, n::Int, P::Tuple{Any, Any})
    old = new = P
    old = (old[1] + 1, old[2])
    while old != new
        old = new
        new = FrobeniusLift(a, h, p, new, n)
    end

    return new
end

function (f::Generic.Frac{<:PolyElem})(x::RingElem)
    return numerator(f)(x)//denominator(f)(x)
end


###############################################################################
#
#   Shifting
#
###############################################################################

#@doc Markdown.doc"""
#    integral(x::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
#> Return the integral of the power series $x$.
#"""
function Generic.integral(x::AbsSeriesElem{T}) where {T <: RingElement}
   xlen = length(x)
   prec = precision(x) + 1
   prec = min(prec, max_precision(parent(x)))
   if xlen == 0
      z = zero(parent(x))
      set_prec!(z, prec)
      return z
   end
   zlen = min(prec, xlen + 1)
   z = parent(x)()
   fit!(z, zlen)
   set_prec!(z, prec)
   z = setcoeff!(z, 0, zero(base_ring(x)))
   for i = 1:xlen
       z = setcoeff!(z, i, coeff(x, i - 1) // base_ring(x)(i))
   end
   set_length!(z, normalise(z, zlen))
   return z
end

#@doc Markdown.doc"""
#    derivative(x::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
#> Return the derivative of the power series $x$.
#"""
function Generic.derivative(x::AbsSeriesElem{T}) where {T <: RingElement}
   xlen = length(x)
   if 1 >= xlen
      z = zero(parent(x))
      set_prec!(z, max(0, precision(x) - 1))
      return z
   end
   z = parent(x)()
   fit!(z, xlen - 1)
   set_prec!(z, precision(x) - 1)
   for i = 1:xlen - 1
      z = setcoeff!(z, i - 1, i * coeff(x, i))
   end
   return z
end


function padic_evaluate(f::Union{SeriesElem, LaurentSeriesFieldElem}, x::RingElem)
    K = base_ring(f)
    ret = zero(K)
    for i in 0:(precision(f) - 1)-valuation(f)
        ret += coeff(f, i+valuation(f))*(K(x)^(i+valuation(f)))
    end
    return ret + O(K, prime(base_ring(f))^precision(f))
end

#function LocalIntegral(F, tP, tQ)
#    f = integral(F)
#    return padic_evaluate(f, tQ) - padic_evaluate(f, tP)
#end

# Return local coordinates on the curve y^a = h(x) around P = (X,Y) up to t-adic precision N.
function LocalCoords(a, h, N, p, P, pts = [])
    if is_in_branch_disk(a, h, P)
        return LocalCoordsB(a, h, N, p, P, pts)
    elseif is_in_inf_disk(a, h, P)
        return LocalCoordsInf(a, h, N, p, P, pts)
    else
        return LocalCoordsNonB(a, h, N, p, P, pts)
    end
end

# Inf
# x coords of all pts should be cong to lambda mod p
# so x coord of P / x coord of pts[i] is a 1-unit hence has an a-th root
# so we normalize all given pts to have x coord same as P
function LocalCoordsInfFancy(a, h, N, p, P, pts = [])
    # below all we really need is the renormalised z-coord but we compute the points anyway for fun.
    npts = [] # normalised points
    for Q in pts
        b = degree(h)
        sc = root(P[1] // Q[1], a) # scale factor
        if length(Q) == 3
            push!(npts, (sc^a * Q[1], sc^b * Q[2], sc * Q[3]))
        else
            # probably the user supplied a point (x,y) so assume z = 1,
            # hope that normalising fixes the valuation
            push!(npts, (sc^a * Q[1], sc^b * Q[2], sc))
        end
    end
    @assert valuation(discriminant(h)) == 0
    K = base_ring(h)
    @assert K.prec_max >= N
    # TODO change to Laurent series once we have derivative
    R,t = LaurentSeriesField(K, N, 't', cached=true)
    #R,t = PowerSeriesRing(K, N, 't', cached=true, model=:capped_absolute)
    # Initial approx
    xt = R(K(P[1]))
    zt = R(K(P[3])) + t
    yt = R(K(P[2]))
    x = gen(parent(h))
    # Newton
    correct_digits = 1
    while correct_digits <= N
        yt = (1//K(a))*(K(a-1)*yt + (reverse(h(P[1]*x))(zt^a))*inv(yt)^(a-1))
        correct_digits *= 2
    end
    if length(pts) > 0
        return (xt,yt,zt, [K(Q[3] - P[3]) for Q in npts])
    end
    return (xt,yt,zt)
end


# This is a different way of doing local coords, we take a lower-brow point of view
# find l,k such that b l - a k = 1 ( as  gcd(a,b) = 1 ).
# Then y has a pole of order b at infinity and x a poly of order a, hence
# t = x^k / y^l is of pole order a k - b l = -1 so is a uniformizer.
# now we solve for x,y using
# y^a - h(x)
# and
# y^l * t - x^k
# simultaneously using Newton, the iteration is
#  [ - h'(x)       | a*y^(a-1)     ] ^ -1
#  [ - k*x^(k-1)   | l*y^(l-1) * t ]
function LocalCoordsInf(a, h, N, p, P, pts = [])
    @assert valuation(discriminant(h)) == 0
    K = base_ring(h)
    b = degree(h)
    _,l,k = gcdx(b,a)
    lambda = lead(h)
    k = -k
    @assert b* l - a* k == 1
    @assert K.prec_max >= N
    # TODO change to Laurent series once we have derivative
    R,t = LaurentSeriesField(K, N, 't', cached=true)
    #R,t = PowerSeriesRing(K, N, 't', cached=true, model=:capped_absolute)
    # Initial approx
    xt = lambda^(-l)*t^(-a) #R(K(P[1]))
    yt = lambda^(-k)*t^(-b) #R(K(P[2]))
    # Newton
    correct_digits = 1

    # TODO maybe redo this with a vector xt,yt, vector f 
    while correct_digits <= N
        M = inv(matrix(R,2,2,[ -derivative(h)(xt), a*yt^(a-1),
                              -k*xt^(k-1), l * yt^(l-1) * t]))

        xtn = xt - (M[1,1]*(yt^a - h(xt)) + M[1,2]*(yt^l * t - xt^k))
        ytn = yt - (M[2,1]*(yt^a - h(xt)) + M[2,2]*(yt^l * t - xt^k))
        xt = xtn
        yt = ytn
        correct_digits *= 2
    end
    @assert xt^k == t * (yt^l)

    if P != :inf
        t0 = K(P[1])^k//K(P[2])^l
    else
        t0 = K(0)
    end
    xt = xt(t)
    yt = yt(t)
    if length(pts) > 0
        return (xt,yt,t0, [K(Q[1])^k//K(Q[2])^l for Q in pts])
    end
    return (xt,yt,t0)
end

# Non-branch
function LocalCoordsNonB(a, h, N, p, P, pts = [])
    @assert valuation(discriminant(h)) == 0
    K = base_ring(h)
    @assert K.prec_max >= N
    # TODO change to Laurent series once we have derivative
    R,t = LaurentSeriesField(K, N, 't', cached=true)
    #R,t = PowerSeriesRing(K, N, 't', cached=true, model=:capped_absolute)
    # Initial approx
    xt = R(P[1]) + t
    yt = R(P[2])
    # Newton
    correct_digits = 1
    while correct_digits <= N
        yt = (1//K(a))*(K(a-1)*yt + (h(xt))*inv(yt)^(a-1))
        correct_digits *= 2
    end
    if length(pts) > 0
        return (xt,yt,K(0), [K(Q[1] - P[1]) for Q in pts])
    end
    return (xt,yt, K(0))
end

# Branch
function LocalCoordsB(a, h, N, p, P, pts = [])
    K = base_ring(h)
    @assert !pos_val(K, discriminant(h))
    #@assert K.prec_max >= N
    # TODO change to Laurent series once we have derivative
    R,t = LaurentSeriesField(K, N, 't', cached=true)
    #R,t = PowerSeriesRing(K, N, 't', cached=true, model=:capped_absolute)
    # Initial approx
    xt = R(P[1])
    yt = R(P[2]) + t
    # Newton
    correct_digits = 1
    while correct_digits <= N
        xtn = xt + (yt^a  - h(xt))*inv(derivative(h)(xt))
        xt = xtn
        correct_digits *= 2
    end
    if length(pts) > 0
        return (xt,yt, K(0), [K(Q[2] - P[2]) for Q in pts])
    end
    return (xt,yt, K(0))
end

function pos_val(K, x)
    return valuation(K(x)) > 0
end

function is_in_branch_disk(a, h, P)
    K = base_ring(h)
    return pos_val(K, P[2])
end

function is_in_inf_disk(a, h, P)
    K = base_ring(h)
    if length(P) == 3
        return pos_val(K, P[3]^degree(h)//P[2])
    end
    return pos_val(K, 1//P[2])
end

function in_same_disk(a, h, P, Q)
    K = base_ring(h)
    return (pos_val(K, P[1] - Q[1]) && pos_val(K, P[2] - Q[2]))
end

function eval_homogeneous(a, h, x, z)
   i = degree(h)
   K = base_ring(h)
   b = degree(h)
   re = zero(K)
   while i >= 0
       re += coeff(h, i)*z^(a*(b-i))*x^(i)
      i -= 1
   end
   return re
end

function verify_pt(a, h, pt)
    if length(pt) == 2
        return pt[2]^a == h(pt[1])
    elseif length(pt) == 3
        b = degree(h)
        return (pt[2])^a == eval_homogeneous(a,h, pt[1], pt[3])
    end
    error("point badly formed")
end

function verify_pts(a, h, pts)
    return all([verify_pt(a, h, p) for p in pts])
end

# returns a p-adic point (X,Y) on y^a = h(x) with x-coord reducing to X,
# and y-coord reducing to y
function lift_point(a, h, X, Y)
    if valuation(Y) > 0
        return lift_y(a, h, Y, X)
    elseif valuation(Y) < 0
        error("not implemented")
    else
        return lift_x(a, h, X, Y)
    end
end

# returns a p-adic point (X,Y) on y^a = h(x) with x-coord X, using hensels lemma
function lift_x(a, h, X, y = nothing)
    K = base_ring(h)
    N = K.prec_max
    if degree(K) == 1
        R = GF(Int(prime(K)))
    else
        R,u = FlintFiniteField(prime(K),degree(K),"u")
    end
    if y == nothing
        y = K(lift_elem(root(R(lift_elem(h(X))), a)))
    end
    if is_in_branch_disk(a,h,(X,y))
        error("not implemented")
    else
        # TODO: in fact this is the same newton iteration as above, can we simplify?
        correct_digits = 1
        while correct_digits <= N
            y = (1//K(a))*(K(a-1)*y + (h(X))*inv(y)^(a-1))
            correct_digits *= 2
        end
    end
    return (X,y)
end

# returns a p-adic point (X,Y) on y^a = h(x) with y-coord Y, using hensels lemma
function lift_y(a, h, Y, x = nothing)
    K = base_ring(h)
    N = K.prec_max
    if x == nothing
        #x = K(lift_elem(root(GF(Int(prime(K)))(lift_elem(h(X))), a)))
        error("not implemented")
    end
    if !is_in_branch_disk(a,h,(x,Y))
        error("not implemented")
    else
        # TODO: in fact this is the same newton iteration as above, can we simplify?
        correct_digits = 1
        while correct_digits <= N
            xn = x + (Y^a  - h(x))*inv(derivative(h)(x))
            x = xn
            correct_digits *= 2
        end
    end
    return (x,Y)

end

# returns a p-adic point (X,Y,Z) on y^a = h(x,z) with z-coord Z, using hensels lemma
function lift_z(a, h, Z, x, y)
    K = base_ring(h)
    N = K.prec_max
    if x == nothing
        #x = K(lift_elem(root(GF(Int(prime(K)))(lift_elem(h(X))), a)))
        error("not implemented")
    end
    varxz = gen(parent(h))
    H = reverse(h(x*varxz))(varxz)
    if !is_in_inf_disk(a,h,(x,y,Z))
        error("not implemented")
    else
        correct_digits = 1
        while correct_digits <= N
            xn = x + (y^a  - H(x))*inv(derivative(H)(x))
            x = xn
            correct_digits *= 2
        end
    end
    return (x,Y)
end

# formal coleman integral \int_P^t in disk around P, if Q supplied return t coord of Q
function FormalTinyColemanIntegralMonomial(a, h, N, p, n, P, i, j; Q = nothing)
    if Q == nothing
        xt,yt,Pt = LocalCoords(a, h, N, p, P)
    else
        xt,yt,Pt,Qt = LocalCoords(a, h, N, p, P, [Q])
    end
    F = xt^i*derivative(xt)*inv(yt)^j
    if Q == nothing
        return integral(F), Pt
    else
        return integral(F), Pt, Qt[1] # only 1 point
    end
end

function TinyColemanIntegralMonomial(a, h, N, p, n, P, Q, i, j)
    @assert in_same_disk(a, h, P, Q)
    K = base_ring(h)
    F,Pt,Qt = FormalTinyColemanIntegralMonomial(a, h, N, p, n, P, i, j; Q = Q)
    return padic_evaluate(F, Qt) - padic_evaluate(F, Pt)
end

function BasisMonomials(a, h)
    return  [(i,j) for j in 1:(a-1) for i in 0:(degree(h)-2)]
end


# Indices of the list from BasisMonomials which are regular 1-forms
# See  Branch Points on Cyclic Covers of the Projective Line - Christopher Towse : Proposition 2.
function RegularIndices(a, h)
    BM = BasisMonomials(a, h)
    k = ceil(degree(h)//a)
    e = k*a - degree(h)
    return [n for n in 1:length(BM) if BM[n][1] < (k*BM[n][2] - 1 - floor(BM[n][2]*e//a))]
end

function FormalTinyColemanIntegralsOnBasis(a::Int, h, N::Int, p::Int, n::Int, P::Tuple)
    return [FormalTinyColemanIntegralMonomial(a, h, N, p, n, P, i, j) for (i,j) in BasisMonomials(a, h)]
end

function TinyColemanIntegralsOnBasis(a::Int, h, N::Int, p::Int, n::Int, P::Tuple, Q::Tuple)
    return elem_type(base_ring(h))[TinyColemanIntegralMonomial(a, h, N, p, n, P, Q, i, j) for (i,j) in BasisMonomials(a, h)]
end

#function rational_points(a, h::Nemo.Poly{Nemo.FinFieldElem})
#    x = FlintQQ(0)
#    ret = Tuple{fmpq,fmpq}[]
#    while height(x) <= bound
#        try
#            y = root(h(x), a)
#            push!(ret, (x,y))
#            if a == 2
#                push!(ret, (x,-y))
#            end
#        catch e
#        end
#        x = next_signed_minimal(x)
#    end
#    return ret
#
#    K = parent(r)
#    u = gen(K)
#    p = characteristic(K)
#    n = degree(K)
#    ui = [u^i for i in 0:(n-1)]
#    for xp in collect(Base.product([1:Int(p) for i in 0:(n-1)]...))
#        x = sum([xp[i] * ui[i] for i in 1:(n-1)])
#        if root(
#            return K(x)
#        end
#    end
#end

function rational_points(a, h::fmpq_poly, bound)
    x = FlintQQ(0)
    ret = Tuple{fmpq,fmpq}[]
    while height(x) <= bound
        try
            y = Coleman.root(h(x), a)
            push!(ret, (x,y))
            if a == 2
                push!(ret, (x,-y))
            end
        catch e
        end
        x = next_signed_minimal(x)
    end
    return ret
end

# TODO wrap proper mult order functions from flint into nemo
# for gfp also nmod?
function multiplicative_order(a)
    N = modulus(parent(a))
    if a == 1
        return 1
    end
    t = a
    i = 1
    while i < N && t != 1
        t *= a
        i+=1
    end
    return i
end

# returns the smallest extension of a padic field containing
# all nth roots of unity.
function Kmun(K, n)
    T = ResidueRing(ZZ, ZZ(n))
    p = prime(K)
    # need n % p^m - 1
    return FlintQadicField(prime(K), multiplicative_order(T(p)), K.prec_max)
end

# vector of Coleman functons I(P) = int_P^y in the disk around x
function ColemanFunctions(a, h, N, p, n, x, y = :inf; frobact = nothing)
    C = ColemanIntegrals(a, h, N, p, n, x, y; frobact = frobact)
    actN = precision(C[1,1])
    F = FormalTinyColemanIntegralsOnBasis(a, h, actN, p, n, x)
    #CastBaseMatrix(parent(C), matrix(base_ring(h), length(tinyints), 1, tinyints))
    K = base_ring(F[1][1])

    return [K(lift_elem(C[i,1])) + O(K, prime(K)^precision(C[i,1])) + F[i][1] for i in 1:length(F)]
end


# Compute coleman ints of a divisor represented as X_1 + X_2 + X_3 ... - Y_1 - Y_2 - Y_3 - ...
# if the second vector Y is unspecified defaults to infinity
function ColemanIntegrals(a, h, N, p, n, x::Array, y = nothing; frobact = nothing)
    if y == nothing
        y = [:inf for P in x]
    end
    return sum([ColemanIntegrals(a, h, N, p, n, x[i], y[i], frobact = frobact) for i in 1:length(x)])
end

function ColemanIntegrals(a, h, N, p, n, y::Symbol, x::Tuple; frobact = nothing)
    return -ColemanIntegrals(a, h, N, p, n, x, y, frobact = frobact)
end

# TODO compute the frobenius action only once
# \int_x^y omega
function ColemanIntegrals(a, h, N, p, n, x::Tuple, y = :inf; frobact = nothing, algorithm = :notdumb)
    if y != :inf
        # Decompose as two integrals, one x to infinity, one y to infinity
        A,B = ColemanIntegrals(a, h, N, p, n, x, :inf) , ColemanIntegrals(a, h, N, p, n, y, :inf)
        #@info A,B
        # original integral is difference of the above
        return A - B
    else
        if is_in_branch_disk(a, h, x)
            # TODO check signs
            if algorithm != :notdumb
                l = Int(divexact(ZZ(p^n - 1),ZZ(a)))
                K = Kmun(base_ring(h), a)
                R,X = PolynomialRing(K, "x")
                D = Hecke.Hensel_factorization(divexact(X^a - 1, X - 1))
                #@info D
                zeta = -coeff([D[k] for k in keys(D) if degree(D[k]) == 1][1],0)
                #zeta = prim_root(K)^(divexact(ZZ(p^degree(K) - 1),a))#teichmuller(K(lift_elem(gen(FiniteField(prime(K), l, "w")[1]))))^l
                @assert zeta^a == 1
                #ex = (p^l - 1)//(p^n - 1) # a gen of F_p^l to the ex power will gen F_p^n
                hK = R([lift_elem(coeff(h,i)) for i in 0:length(h)])
                # TODO check this is always compatible, not always conway!! for large p
                xK = (K(lift_elem(x[1])),K(lift_elem(x[2]))) # x base changed to K
                muxK = superelliptic_automorphism(a, hK, p, degree(K), xK)
                #@info xK
                #@info muxK

                return [get_padic(inv(zeta^(-ij[2]) - 1)*b) + O(base_ring(h), prime(base_ring(h))^N) for (ij, b) in zip(BasisMonomials(a, h), TinyColemanIntegralsOnBasis(a, hK, N, p, degree(K), xK, muxK))]
            else
                B = lift_y(a, h, base_ring(h)(0), x[1])
                #@info x
                #@info B

                return TinyColemanIntegralsOnBasis(a, h, N, p, n, B, x)
            end
        else
            pts = [x]
            for i in 1:n-1
                push!(pts, FrobeniusLift(a, h, p, pts[end]))
            end
            #@info "pts",pts
            M, C = AbsoluteFrobeniusActionOnLift(a, h, N, p, n, pts)
            #@info "bigC",C
            #@info("tt")
            #@info M, C
            # C is the matrix ( f_i(P) | f_i(phi(P)) | --- | f_i(phi^{n-1}(P)) )

            # get the q-power action see Remark 12 of Balakrishnan-Bradshaw-Kedlaya
            prodphiM = identity_matrix(M)
            CC = zero_matrix(base_ring(C), nrows(C), 1)
            #@info "M",M
            #@info "CC",CC
            frobMs = [M]
            for i in 1:n-1
                push!(frobMs, map(frobenius,frobMs[end]))
            end
            #@info "frobM",frobMs
            for t in n-1:-1:0
                #@info "t",t
                #@info "CC",CC
                #@info "M",prodphiM
                #@info "Ccol",C[1:nrows(C), t + 1]
                #@info "fM",map(frobenius,prodphiM)
                #@info "fCcol",map(frobenius, C[1:nrows(C), t + 1])
                CC += prodphiM * map(x->frobenius(x, t), C[1:nrows(C), t + 1])#0 ? frobenius : x->x, C[1:nrows(C), t + 1])
                prodphiM = prodphiM * frobMs[t + 1]
            end


            #@info(x)
            #@info(FrobeniusLift(a, h, p, x, n))
            tinyints = TinyColemanIntegralsOnBasis(a, h, N, p, n, x, FrobeniusLift(a, h, p, x, n))
            #@info("ti: \n",tinyints)
            #@info("C: \n",C)
            #@info("M\n", trace(M))
            #@info("M-1\n", M-1)
            #@info("inv(M-1)\n", inv(M-1))
            #@info("inv(M-1)C\n", inv(M-1)*C)
            #@info("M-1\n", M-1)
            #@info("diff:\n",C - CastBaseMatrix(parent(C), matrix(base_ring(h), length(tinyints), 1, tinyints)))

            #@info tinyints
            #@info inv(prodphiM - 1) * (CC - CastBaseMatrix(parent(CC), matrix(base_ring(h), length(tinyints), 1, tinyints)))
            return inv(prodphiM - 1) * (CC - CastBaseMatrix(parent(CC), matrix(base_ring(h), length(tinyints), 1, tinyints)))
        end
    end
end

function weierstrassprep(f::Generic.RelSeries)
    lowval = 2<<30
    index = -1

    for i in valuation(f):precision(f)
        if valuation(coeff(f, i)) < lowval
            lowval = valuation(coeff(f, i))
            index = i
        end
    end
    @assert lowval == 0
    R,t = PolynomialRing(base_ring(f), 't')
    g = zero(R)
    g2 = zero(parent(f))
    o = f
    for i in 1:base_ring(f).prec_max
        o = o//shift_right(o, index)
    end

    return (o, f//o)
end

function residue_field(R::FlintPadicField)
    return FlintFiniteField(prime(R), 1, "w")[1]
end

function residue_field(R::FlintQadicField)
    return FlintFiniteField(prime(R), degree(R), "w")[1]
end

function prim_root(K::FlintPadicField)
    return teichmuller(K(coeff(lift_elem(gen(residue_field(K))),0)))
end

function prim_root(K::FlintQadicField)
    return teichmuller(K(lift_elem(gen(residue_field(K)))))
end

# how many roots (via strassman)
function strassman_num_roots(f)
    p = prime(base_ring(f))
    @assert p != 2
    curmin = valuation(coeff(f, 0))
    n = valuation(f)
    N = valuation(f)
    while n < curmin*(p-1)//(p-2)
        if valuation(coeff(f, n) + n) <= curmin
            if n < precision(f)
                curmin = valuation(coeff(f, n) + n)
                N = n
            else
                return 10000
            end
        end
        n += 1
    end
    return N
end

# TODO compute frob as little as possible
function Chabauty(a, h, N, p, n, jacobian_gens, res_disks = nothing)
    K = base_ring(h)
    if res_disks == nothing
        kappa = residue_field(K)
        Rb,x = PolynomialRing(kappa, "x")
        hb = Rb(lift_elem(h))
        res_disks = rational_points(a, hb)
    end
    r = length(jacobian_gens)
    g = divexact(((a-1)*(degree(h)-1)),2)
    Cols = [ColemanIntegrals(a, h, N, p, n, P, Q) for (P,Q) in jacobian_gens]
    K2 = K
    if length(Cols) > 0
        K2 = parent(Cols[1][1,1]) # probably this is K, but maybe diff prec
    end
    RI = RegularIndices(a, h)
    #@info(Cols)
    colsmat = matrix(K2,[Cols[j][i,1] for j in 1:length(jacobian_gens), i in RI])
    #@info(colsmat)
    _,annih = nullspace(colsmat)
    #@info("aani",annih)
    pts = []
    for Pb in res_disks
        P = lift_point(a, h, K(lift_elem(Pb[1])), K(lift_elem(Pb[2])))
        @assert verify_pt(a, h, P)
        cfs = ColemanFunctions(a, h, N, p, n, P, :inf) # TODO this could be made faster by passing RI in
        regcfs = [cfs[i] for i in RI]
        Kf = base_ring(regcfs[1])
        local_ann_funcs = [sum([(Kf(lift_elem(annih[i,j])) + O(Kf, p^precision(annih[i,j])))*regcfs[i] for i in 1:g]) for j in 1:ncols(annih)]
        if sum([sum([valuation(coeff(l,i)) for i in 0:precision(l)-1]) for l in local_ann_funcs]) > 0 || true
            #@info(P)
            #@info("lafs ", local_ann_funcs)
            #@info("lafs ", [[valuation(coeff(l,i)) for i in 0:precision(l)-1] for l in local_ann_funcs])
            #@info("nr ", [strassman_num_roots(f) for f in local_ann_funcs])
            #@info("\n\n")
        end
    end
    return pts
end

end # module
