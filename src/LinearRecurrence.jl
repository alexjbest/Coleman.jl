#  Copyright (C) 2009-2011 Moriz Minzlaff <minzlaff@daad-alumni.de>
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
# This package implements Bostan et al.'s algorithm for linear
# recurrences [1] with modifications a la Harvey [2].
#
# [1] Bostan, A., Gaudry, P., Schost, E.: Linear recurrences with polynomial
#   coefficients and application to integer factorization and Cartier-Manin
#   operator. SIAM J. Comput. 36, 1777–1806 (2007)
#
# [2] Harvey, D.: Kedlaya's algorithm in larger characteristic. Int. Math.
#   Res. Notices 2007 (2007)
#
# [3] von zur Gathen, J., Gerhard J.: Modern Computer Algebra. Cambridge
#   University Press (2003)
#
#*****************************************************************************

using Nemo

import AbstractAlgebra.PolyElem
import AbstractAlgebra.SeriesElem

function xadic_to_padic(f::SeriesElem, R::PadicField)
    ret = zero(R)
    p = prime(R)
    for i in 0:(precision(f) - 1)-valuation(f)
        ret += coeff(f, i)*p^(i+valuation(f))
    end
    return ret + O(R,prime(R)^precision(f))
end

function (R::FlintIntegerRing)(n::fmpz_poly)
    return R(coeff(n,0))
end

#function (R::FlintIntegerRing)(n::fmpq)
#    @assert denominator(n) == 1
#    return R(numerator(n))
#end

function (R::FlintIntegerRing)(n::fmpq_poly)
    return R(coeff(n,0))
end

function cast_poly_nmod(R, e)
    RR = base_ring(R)
    # If K/L with L prime, then RR is the UnramifiedQuotientRing(L,N)
    # If R is not a polynomial ring, then (RR eq RRR)
    res_ = Vector{elem_type(RR)}(undef, degree(e) + 1)#[zero(RR) for i in 1:(degree(e) + 1)]
    for i = 0:degree(e)
        res_[i+1] = RR(lift_elem(coeff(e,i)))
    end
    return R(res_)
end

function lift_elem(ei)
    return lift(FmpzPolyRing(FlintZZ,:x), ei)
end

function lift_elem(ei::Integer)
    return ei
end

function lift_elem!(r::fmpz_poly, a::qadic)
  ctx = parent(a)
  res = Bool(ccall((:padic_poly_get_fmpz_poly, :libflint), Cint,
                   (Ref{fmpz_poly}, Ref{qadic}, Ref{FlintQadicField}), r, a, ctx))
  !res && error("Unable to lift")
  return r
end

function lift_elem(ei::Nemo.FinFieldElem)
    R = FmpzPolyRing(FlintZZ,:x)
    x = gen(R)
    o = zero(R)
    for i in 0:degree(parent(ei))-1
        o += coeff(ei, i)*x^i
    end
    return o
end

function lift_elem(ei::Nemo.padic)
    return lift(FlintQQ, ei)
end

function lift_elem!(r::fmpq, a::padic)
  ctx = parent(a)
  ccall((:padic_get_fmpq, :libflint), Nothing,
        (Ref{fmpq}, Ref{padic}, Ref{FlintPadicField}), r, a, ctx)
  return r
end

function lift_elem(ei::Nemo.fmpq)
    return ei
end

function lift_elem(ei::Nemo.gfp_elem)
    return ei.data
end

function lift_elem(ei::Generic.ResElem{<:Union{RingElem, Integer}})
    return data(ei)
end

function LowerCaseDD(alpha,beta,d)
    #    Return the product dd(alpha,beta,d)
    #
    #    INPUT:
    #
    #    -  ``alpha`` - an element of R
    #    -  ``beta`` - an element of R
    #    -  ``d`` - a finite cardinal
    #
    #    OUTPUT:
    #
    #    Returns the element dd(alpha,beta,d) of R
    #
    res = beta
    #for i = 2:d
    #    j = j * i
    #    #res = res * i
    #    @assert 
    #end
    res = beta * factorial(fmpz(d))
    for i = -d:d
      t = parent(beta)(i)
      t = mul!(t, t, beta)
      t = add!(t, alpha, t)
      res = mul!(res, res, t)
    end
    return res
end

function Nemo.add!(z::fmpq_abs_series, a::fmpq_abs_series, b::fmpq_abs_series)
   lena = length(a)
   lenb = length(b)

   prec = min(a.prec, b.prec)

   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   ccall((:fmpq_poly_add_series, :libflint), Nothing,
                (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Ref{fmpq_abs_series}, Int),
               z, a, b, lenz)
   return z
end


function LowerCaseDD_(alpha, beta, d)
    #    Return all factors of dd(alpha,beta,d)
    #
    #    INPUT:
    #
    #    -  ``alpha`` - an element of R
    #    -  ``beta`` - an element of R
    #    -  ``d`` - a finite cardinal
    #
    #    OUTPUT:
    #
    #    Returns a sequence whose entries are the factors of dd(alpha,beta,d)
    #    ordered as before Theorem 5
    #
    R = parent(beta)
    res_ = Vector{elem_type(R)}(undef, 3*d + 1)#[parent(beta)(0) for i in 1:(3*d + 1)]
    res_[1] = beta
    for i = 2:d
        res_[i] = R(i)
    end

    if d > 1
      res_[2 * d + 1 - d] = alpha - res_[d] * beta
    else
      res_[2 * d + 1 - d] = alpha - R(1) * beta
    end
    for i = (-d + 1):d
      res_[2*d+1+i] = res_[2*d + 1 + (i - 1)] + beta# i*beta
    end

    return res_
end



function UpperCaseDD_(alpha, beta, k)
    #    Return all factors of D(alpha,beta,k)
    #
    #    INPUT:
    #
    #    -  ``alpha`` - an element of R
    #    -  ``beta`` - an element of R
    #    -  ``k`` - a finite cardinal
    #
    #    OUTPUT:
    #
    #    Returns a sequence whose entries are the factors of D(alpha, beta, k) as
    #    defined before Lemma 6
    #
    k_ = Int[ k ]
    while k_[end] > 1
        push!(k_, (k_[end] >> 1))
    end

    R = parent(alpha)

    res_ = Vector{elem_type(R)}(undef, 2 * (length(k_) - 1))
    for i = 1:(length(k_)-1)
        res_[2*i-1] = LowerCaseDD(beta*(k_[i+1]+1),beta,k_[i+1])
        res_[2*i] = LowerCaseDD(alpha*k_[i+1],beta,k_[i+1])
    end
    return res_
end

function RetrieveInverses(prodInv,r_)
    #    Implements [1, Lemma 1]
    #
    #    INPUT:
    #
    #    -  ``prodInv`` - the inverse of the product of the elements in r_
    #    -  ``r_`` - sequence of ring elements
    #
    #    OUTPUT:
    #
    #    A sequence inv_ such that inv_[i] = r_[i]^{-1}
    #
    if (length(r_) == 0)
        return typeof(prodInv)[prodInv]
    end
    R = parent(r_[1])
    #rProd = [parent(r_[1])(0) for i in 1:length(r_)]
    rProd = Vector{elem_type(R)}(undef, length(r_))
    rProd[1] = r_[1]
    for i = 2:length(r_)
        rProd[i] = r_[i]*rProd[i-1]
    end
    #inv_ = [parent(prodInv)(0) for i in 1:length(r_)]
    inv_ = Vector{elem_type(R)}(undef, length(r_))
    inv_[length(r_)] = deepcopy(prodInv)
    for i = (length(r_)-1):-1:1
        inv_[i] = r_[i+1]*inv_[i+1]
    end
    for i = 2:length(r_)
        inv_[i] = mul!(inv_[i], inv_[i], rProd[i-1])
    end
    return inv_
end

function better_evaluate(a::PolyElem{T}, b::T) where {T}
  i = length(a)
  R = parent(b)
  if i == 0
    return zero(R)
  end
  z = one(R)
  z = mul!(z, z, coeff(a, i - 1))
  while i > 1
    i -= 1
    z = mul!(z, z, b)
    add!(z, z, coeff(a, i - 1))
  end
  return z
end

function better_evaluate(a, b)
  return evaluate(a, b)
end

function Evaluate(M, x)
    return map(t->better_evaluate(t, x), M)
end

function ShiftEvaluationPre(alpha, beta, ddi_, d, RPol)
    #    Implements precompuations for [1, Theorem 5]
    #
    #    INPUT:
    #
    #    -  ``alpha`` - an element of R
    #    -  ``beta`` - an element of R
    #    -  ``ddi_`` - contains the inverses of the factors of
    #                 dd(alpha, beta, d) in the order of the definition
    #    -  ``d`` - a finite cardinal (= #baseValues_)
    #    -  ``RPol`` - univariate polynomial ring over R
    #
    #    OUTPUT:
    #
    #    Returns the two sequences partiali_ and delta_ such that
    #    partiali_[i] = partial(i-1,d)^{-1} and delta_[i] =
    #    delta(a, i-1, d) as in BGS07, Lemma 2
    #
    R = base_ring(RPol)
    #partiali_ = [zero(R) for i in 1:(d + 1)]
    partiali_ = Vector{elem_type(R)}(undef, d + 1)
    partiali_[1] = one(R)
    for i = 2:d
        partiali_[1] = mul!(partiali_[1], partiali_[1],  ddi_[i])
    end
    if (ModByPowerOf2(d,1) == 1)
        partiali_[1] = -partiali_[1]
    end
    partiali_[2] = -d*partiali_[1]
    for i = 3:d+1
        partiali_[i] = R((i-1)-d-1)
        partiali_[i] = mul!(partiali_[i], partiali_[i], ddi_[i-1])
        partiali_[i] = mul!(partiali_[i], partiali_[i], partiali_[i-1])
    end

    a = alpha*ddi_[1]
    delta_ = Vector{elem_type(R)}(undef, d + 1)
    delta_[1] =  one(R)
    for i = 0:d
        delta_[1] = mul!(delta_[1], delta_[1], a - i)
    end
    for i = 2:d+1
      delta_[i] = (a+i-1)*(ddi_[d+i-1])
      delta_[i] = mul!(delta_[i], delta_[i], beta)
      delta_[i] = mul!(delta_[i], delta_[i], delta_[i-1])
    end

    s = RPol(elem_type(R)[ ddi_[d+i]*beta for i in 1:(2*d+1) ])

    return partiali_, delta_, s
end

function ShiftEvaluation(partiali_, delta_, s, ddi_, d,
                         baseValues_, RPol)
    #    Implements [1, Theorem 5]
    #
    #    INPUT:
    #
    #    -  ``partiali_`` - a sequence of elements of a ring R
    #    -  ``delta_`` - a sequence of elements of R
    #    -  ``s`` - a polynomial over R
    #    -  ``ddi_`` - contains the inverses of the factors of
    #                 dd(alpha, beta, d) in the order of the definition in BGS07
    #    -  ``d`` - a finite cardinal
    #    -  ``baseValues_`` - a sequence with i-th entry is f(r + (i-1)*beta)
    #                         where f is a polynomial over R of degree d
    #    -  ``RPol`` - polynomial ring over R
    #
    #    OUTPUT:
    #
    #    Returns a sequence with i-th entry f(r + (i-1)*beta + alpha)
    #
    #    NOTE:
    #    Complexity: 2M(d) + O(d) = O( M(d) )
    #    (together with ShiftEvaluationPre)
    #
    R = base_ring(RPol)
    p = RPol(elem_type(R)[ baseValues_[i]*partiali_[i] for i in 1:(d+1) ])
    #q = p*s   # this multiplication accounts for roughly
    # 1/3 of all computation time spent in LinearRecurrence
    #res_ = elem_type(R)[ R(delta_[k+1]*coeff(q,d+k)) for k in 0:d ]
    res_ = _coeffs_in_the_middle(p, s, d)
    for i in 1:(d + 1)
      res_[i] = mul!(res_[i], res_[i], delta_[i])
    end

    return res_
end

# Given f of degree d and g of degree 2*d, return
# [ coeff(f * g, d + k) for k in 0:d ]
# without computing f * g
function _coeffs_in_the_middle(f, g, d)
  R = base_ring(f)
  t = R()
  z = elem_type(R)[zero(R) for i in 0:d]
  for k in d:(2*d)
    for j in 0:k 
      t = mul!(t, coeff(g, j), coeff(f, k - j))
      z[k - d + 1] = add!(z[k - d + 1], z[k - d + 1], t)
    end
  end
  return z
end

function Algorithm10_3(moduli_::Vector{T}, k) where {T}
    #    Implements [3, Algorithm 10.3]
    #
    #    NOTE: #moduli_ = 2^k
    #
    cache = Vector{T}[ T[ zero(parent(moduli_[1])) for i in 1:2^(j - 1) ] for j in (k + 1):-1:1]
    n = length(moduli_)
    res_ = Vector{T}[T[] for i in 1:k+1]
    res_[1] = T[ moduli_[j] for j in 1:n]
    for j in 1:n
      cache[1][j] = set!(cache[1][j], moduli_[j])
    end
    for i = 2:k+1
        for j = 1:(n >> (i-1))
            push!(res_[i], res_[i-1][2*j-1]*res_[i-1][2*j])
            cache[i][j] = mul!(cache[i][j], res_[i - 1][2 * j - 1], res_[i - 1][2 * j])
        end
    end
    #@show length.(res_)
    #@show length.(cache)
    @assert cache == res_
    return res_
end

function Algorithm10_3_cached(moduli_::Vector{T}, k, cache) where {T}
    #    Implements [3, Algorithm 10.3]
    #
    #    NOTE: #moduli_ = 2^k
    #
    n = length(moduli_)
    #res_ = Vector{T}[T[] for i in 1:k+1]
    #res_[1] = T[ moduli_[j] for j in 1:n]
    for j in 1:n
      cache[1][j] = set!(cache[1][j], moduli_[j])
    end
    for i = 2:k+1
        for j = 1:(n >> (i-1))
            #push!(res_[i], res_[i-1][2*j-1]*res_[i-1][2*j])
            #cache[i][j] = mul!(cache[i][j], cache[i - 1][2 * j - 1], cache[i - 1][2 * j])
            cache[i][j] = cache[i - 1][2 * j - 1] * cache[i - 1][2 * j]
        end
    end
    #@show length.(res_)
    #@show length.(cache)
    #@assert cache == res_
    return cache
end

function Algorithm10_5(f, moduli_, Mij_)
    #    Implements [3, Algorithm 10.5]
    #
    #    NOTE:
    #    #moduli_ = 2^k
    #    #Mij_ = k+1
    #
    n = length(moduli_)
    k = length(Mij_) -1
    #new_rem_ = [ 0 * f for i in 1:n]
    #new_rem2_ = [ 0 * f for i in 1:n]
    #new_rem_[1] = deepcopy(f)

    rem_ = typeof(f)[ f ]
    while (k > 0)
        remOld_ = rem_
        rem_ = Vector{typeof(f)}(undef, 2*(n >> k))
        for i = 1:(n >> k)
            rem_[2*i - 1] = mod(remOld_[i], Mij_[k][2*i-1])
            rem_[2*i]     = mod(remOld_[i], Mij_[k][2*i])
    #        mod!(new_rem2_[2*i - 1], new_rem_[i], Mij_[k][2*i-1])
    #        mod!(new_rem2_[2*i],  new_rem_[i], Mij_[k][2*i])
    #        @assert rem_[2 * i - 1] == new_rem2_[2 * i - 1]
    #        @assert rem_[2 * i] == new_rem2_[2 * i]
        end
        k = k-1
    #    new_rem_, new_rem2_ = new_rem2_, new_rem_
    end
    #if isodd(k)
    #  new_rem_, new_rem2_ = new_rem2_, new_rem_
    #end
    #@assert new_rem_ == rem_

    res_ = elem_type(base_ring(parent(f)))[ coeff(rem_[i],0) for i in 1:length(rem_) ]
    #res2_ = elem_type(base_ring(parent(f)))[ coeff(new_rem_[i],0) for i in 1:length(new_rem_) ]
    #@assert res_ == res2_
    return res_
end

function Algorithm10_5_cached(f, moduli_, Mij_, cache1, cache2)
    #    Implements [3, Algorithm 10.5]
    #
    #    NOTE:
    #    #moduli_ = 2^k
    #    #Mij_ = k+1
    #
    n = length(moduli_)
    k = length(Mij_) -1
    #@assert n == length(cache1)
    #@assert n == length(cache2)
    #cache1[1] = deepcopy(f)
    set!(cache1[1], f)

    #rem_ = typeof(f)[ deepcopy(f) ]
    while (k > 0)
    #    remOld_ = rem_
    #    rem_ = Vector{typeof(f)}(undef, 2*(n >> k))
        for i = 1:(n >> k)
    #        rem_[2*i - 1] = mod(remOld_[i], Mij_[k][2*i-1])
    #        rem_[2*i]     = mod(remOld_[i], Mij_[k][2*i])
            #cache2[2*i - 1] = mod(cache1[i], Mij_[k][2*i-1])
    #        z = zero(parent(f))
    #        mod!(z, cache1[i], Mij_[k][2* i - 1])
    #        @assert z == mod(cache1[i], Mij_[k][2*i-1])
            cache2[2 * i] =mod!(cache2[2*i], cache1[i], Mij_[k][2*i])
            #z = zero(parent(f))
            #mod!(z, cache1[i], Mij_[k][2* i])
            #@assert z == mod(cache1[i], Mij_[k][2*i])
            cache2[2* i - 1] = mod!(cache2[2*i - 1], cache1[i], Mij_[k][2*i-1])
            #@assert cache2[2 * i - 1] == mod(cache1[i], Mij_[k][2*i-1])
            #@assert cache2[2 * i] == mod(cache1[i], Mij_[k][2*i])
            #@assert rem_[2 * i - 1] == cache2[2 * i - 1]
            #@assert rem_[2 * i] == cache2[2 * i]
        end
        k = k-1
        cache1, cache2 = cache2, cache1
    end
    if isodd(k)
      cache1, cache2 = cache2, cache1
    end

    #@assert cache1 == rem_

    res_ = elem_type(base_ring(parent(f)))[ coeff(cache1[i],0) for i in 1:length(cache1) ]
    return res_
end

function Algorithm10_9(points_, cvalues_, Mij_, RPol)
    #    Implements [3, Algorithm 10.9]
    #
    #    NOTE:
    #    #points_ = #values_ = 2^k
    #    #Mij_ = k+1
    #
    k = length(Mij_) -1
    if (k == 0)
        return RPol(cvalues_[1])
    end

    n = length(points_)
    Mij1_ = eltype(Mij_)[ Mij_[i][1:(n >> i)] for i in 1:k ]
    Mij2_ = eltype(Mij_)[ Mij_[i][(n >> i)+1:(n >> (i-1))] for
             i in 1:k ]
    res_ = Mij_[k][2]*Algorithm10_9(points_[1:(n >> 1)],
                                   cvalues_[1:(n >> 1)], Mij1_, RPol) +
            Mij_[k][1]*Algorithm10_9(points_[(n >> 1)+1:n],
                                    cvalues_[(n >> 1)+1:n], Mij2_, RPol)
    return res_
end

function MatrixAPEvaluationPre(k, alpha::T, beta, DDi_, RPol) where {T}
    #    Implements precomputations for [1, Theorem 8]
    #
    #    INPUT:
    #
    #    -  ``k`` - a finite cardinal
    #    -  ``alpha`` - an element of R
    #    -  ``beta`` - an element of R
    #    -  ``DDi_`` - a sequence containing the inverses of the lower case
    #                  dds contained in the upper case DD, the order
    #                  is the same as in the definition before Lemma 6,
    #                  i.e. DDi_[2*i-1] = dd(beta*(k_{i}+1),beta,k_{i})^(-1),
    #                       DDi_[2*i] = dd(alpha*k_{i},beta,k_{i})^(-1)
    #    - ``RPol`` - a polynomial ring over R
    #
    #    OUTPUT:
    #
    #    Returns a sequence of integers k_, an integer logk, and sequences
    #    of sequences of elements of R ddi1__, partiali1__, delta1__,
    #    and a sequence of polynomials s1_ over R (and the same with
    #    2 instead of 1)
    #
    k_ = Int[ k ]             			# k_[i] = k_{i-1}
    while k_[end] > 1              # k_[1] = k_0 = k
        push!(k_, k_[end] >> 1)
    end
    logk = length(k_)-1				# = Floor(Log(2,k))

    ddi1__ = Vector{Vector{T}}(undef, logk)#[typeof(alpha)[] for i in 1:logk]
    partiali1__ = Vector{Vector{T}}(undef, logk)#[typeof(alpha)[] for i in 1:logk]
    delta1__ = Vector{Vector{T}}(undef, logk)#[typeof(alpha)[] for i in 1:logk]
    s1_ = Vector{elem_type(RPol)}(undef, logk)
    ddi2__ = Vector{Vector{T}}(undef, logk)
    partiali2__ = Vector{Vector{T}}(undef, logk)
    delta2__ = Vector{Vector{T}}(undef, logk)
    s2_ = Vector{elem_type(RPol)}(undef, logk)
    for i = logk:-1:1
        d = k_[i+1]
        tbeta = (d + 1) * beta
        talpha = d * alpha
        dd1_ = LowerCaseDD_(tbeta, beta, d)
        ddi1__[i] = RetrieveInverses(DDi_[2*i-1], dd1_)
        dd2_ = LowerCaseDD_(talpha, beta, d)
        ddi2__[i] = RetrieveInverses(DDi_[2*i], dd2_)
        partiali1__[i], delta1__[i], s1_[i] =
        ShiftEvaluationPre(tbeta, beta, ddi1__[i],
                           k_[i+1], RPol)
        partiali2__[i], delta2__[i], s2_[i] =
        ShiftEvaluationPre(talpha, beta, ddi2__[i],
                           k_[i+1], RPol)
    end

    return k_, logk, ddi1__, partiali1__, delta1__, s1_, ddi2__,partiali2__, delta2__, s2_
end

function ModByPowerOf2(n,e)
    return n - (n >> e) << e
end


function MatrixAPEvaluation(M, k_, logk,
                            ddi1__, partiali1__, delta1__, s1_,
                            ddi2__, partiali2__, delta2__, s2_, alpha, beta, DDi_)
    #    Implements [1, Theorem 8]
    #
    #    INPUT:
    #
    #    -  ``M`` - matrix of linear polynomials
    #    -  ``k_`` - a sequence of finite cardinals
    #    -  ``logk`` - a finite cardinal
    #    -  ``ddi1__`` - a sequence of sequences of elements of a ring R
    #    -  ``partiali1__`` - a sequence of sequences of elements of a ring R
    #    -  ``delta1__`` - a sequence of sequences of elements of a ring R
    #    -  ``s1_`` - a sequence of polynomials over R
    #    -  ``ddi2__`` - a sequence of sequences of elements of a ring R
    #    -  ``partiali2__`` - a sequence of sequences of elements of a ring R
    #    -  ``delta2__`` - a sequence of sequences of elements of a ring R
    #    -  ``s2_`` - a sequence of polynomials over R
    #    -  ``alpha`` - an element of R
    #    -  ``beta`` - an element of R
    #    -  ``DDi_`` - a sequence containing the inverses of the lower case
    #                  dds contained in the upper case DD, the order
    #                  is the same as in the definition before Lemma 6,
    #                  i.e. DDi_[2*i-1] = dd(beta*(k_{i}+1),beta,k_{i})^(-1),
    #                       DDi_[2*i] = dd(alpha*k_{i},beta,k_{i})^(-1)
    #    OUTPUT:
    #
    #    Returns a sequence with i-th entry M_k((i-1)*beta) where
    #    M_k(X) = M(X + alpha)...M(X + alpha*k)
    #
    #    NOTE:
    #    Complexity (together with MatrixAPEvaluationPre):
    #    If M is an nxn matrix, then the runtime is O( MM(n)k + n^2M(k) )
    #
    # PRELIMINARIES
    R = base_ring(base_ring(M))
    RPol = base_ring(M)
    x = gen(RPol)
    n = nrows(M)
    alpha = R(alpha)
    beta = R(beta)

    # STEP logk+1
    res_ = dense_matrix_type(R)[ Evaluate(M, alpha), Evaluate(M, alpha+beta) ]
    old1_ = dense_matrix_type(R)[ zero_matrix(R, n, n), zero_matrix(R, n, n) ]
    old2_ = dense_matrix_type(R)[ zero_matrix(R, n, n), zero_matrix(R, n, n) ]

    # STEPS logk:2
    for i = logk:-1:1
        # at the start of step i
        # res_[l] = M_k_[i+1]((l-1)*beta), l = 1,...,k_[i+1]+1
        #
        # during step i we set
        # 	old1_[l] = M_k_[i]((l-1)*beta)
        # 	old2_[l] = M_k_[i]((l-1)*beta + alpha*k_[i])
        #
        # at the end of the step i
        # 	res_[l] = M_k_[i]((l-1)*beta), l = 1,...,k_[i]+1
        #
        old1_ = append!(old1_, (zero_matrix(R, n, n) for l in (length(old1_)+1):(2*(k_[i+1]+1))))
        # make space for new matrices
        old2_ = append!(old2_, ( zero_matrix(R, n, n) for l in (length(old2_)+1):(2*(k_[i+1]+1)) ))

        if (ModByPowerOf2(k_[i],1) == 1)
            # "correction" term in case of odd k_[i]
            correction_ = dense_matrix_type(R)[ Evaluate(M, l*beta + alpha*k_[i]) for
                           l in 0:(k_[i]) ]
        end
        d = k_[i+1]

        baseValues_ = Vector{elem_type(R)}(undef, k_[i+1]+1)

        for r = 1:n
            # to deduce the components of M_k_[i]
            for c = 1:n
                # we need more values of each component of M_k_[i+1]
                for j in 1:(k_[i+1]+1)
                  baseValues_[j] = res_[j][r, c]
                end
                #baseValues_ = elem_type(R)[ res_[j][r,c] for j in 1:(k_[i+1]+1) ]
                #@assert baseValues2_ == baseValues_
                values1_ = ShiftEvaluation(partiali1__[i], delta1__[i], s1_[i], ddi1__[i], d, baseValues_, RPol)
                shiftedValues_ = ShiftEvaluation(partiali2__[i], delta2__[i], s2_[i], ddi2__[i], d, baseValues_, RPol)
                values2_ = ShiftEvaluation(partiali1__[i], delta1__[i], s1_[i], ddi1__[i], d, shiftedValues_, RPol)
                for l in 1:length(baseValues_)
                  old1_[l][r,c] = baseValues_[l]
                  old2_[l][r,c] = shiftedValues_[l]
                end
                kk = length(baseValues_)
                for l in (length(baseValues_)+ 1):2*(k_[i+1]+1)
                  old1_[l][r,c] = values1_[l - kk]
                  old2_[l][r,c] = values2_[l - kk]
                end
            end
        end
        # using the evaluated matrices, compute the new res_
        if (ModByPowerOf2(k_[i],1) == 0)
            res_ = dense_matrix_type(R)[ old1_[l]*old2_[l] for l in 1:(k_[i]+1) ]
        else
            res_ = dense_matrix_type(R)[ old1_[l]*old2_[l]*correction_[l] for
                    l in 1:(k_[i]+1) ]
        end
    end

    return res_
end

function ToBase(M)
    Mb = matrix(base_ring(base_ring(M)),
                elem_type(base_ring(base_ring(M)))[ coeff(M[i, j], 0) for i = 1:nrows(M),j = 1:ncols(M)])
    return Mb
end

function MatrixEvaluationPre(M, k, ki_)
    #    Implements precomputations for [1, Corollary 10]
    #
    #    INPUT:
    #
    #    - ``M`` - a matrix of linear polynomials
    #    - ``k`` - an integer > 0
    #    - ``ki_`` - a sequence of integers such that ki_[i] = i^{-1},
    #                i=1, ..., k2-1 (k2 is the smallest power of 2 >k)
    #
    #    OUTPUT:
    #
    k2 = length(ki_) +1

    # 1: Proposition 9 of BGS07, Evaluate Mk at 0, 1, ..., k2-1
    R = base_ring(base_ring(M))
    n = nrows(M)
    M_ = dense_matrix_type(R)[ Evaluate(M,R(i)) for i in 1:(k+k2-1) ]
    L_ = Vector{dense_matrix_type(R)}(undef, k2)
    for i = k2:-1:(k+1)
        L_[i] = identity_matrix(R, n)
    end
    for i = k:-1:1
        L_[i] = M_[i]*L_[i+1]
    end		# L_ is ok, collects M(1),...,M(k)
    C_ = Vector{dense_matrix_type(R)}(undef, k2)
    C_[1] = identity_matrix(R, n)
    for i = 2:(k+1)
        C_[i] = C_[i-1]*M_[k+i-1]
    end
    C_[k2] = M_[2*k]
    for j = (2*k-1):-1:k2
        C_[k2] = M_[j]*C_[k2]
    end
    for i = (k2-1):-1:(k+2)
        C_[i] = M_[i]*C_[i+1]
    end		# C_ is ok, collects M(k+1),...,M(2*k)
    R_ = Vector{dense_matrix_type(R)}(undef, k2)
    for i in 1:k+1
      R_[i] = identity_matrix(R, n)
    end
    for i = (k+2):k2
        R_[i] = R_[i-1]*M_[k+i-1]
    end		# R_ is ok, collects M(2*k+1),...,M(2*k+(k2-k))

    MkEval_ = Vector{dense_matrix_type(R)}(undef, k2)
    for i in 1:k2
      MkEval_[i] = L_[i]*C_[i]
      MkEval_[i] = mul!(MkEval_[i], MkEval_[i], R_[i])
    end
    # at  this point, MkEval_[i] contains M_k(i-1), 1 \le i \le k2

    # 2: Interpolate:compute M_k(X) a la Algorithm 10.11
    # zur Gathen/Gerhard
    # i. this replaces the use of algorithm 10_5
    facki_ = Vector{elem_type(R)}(undef, k2 - 1) # facki_[i] contains 1/fac(i)
    facki_[1] = one(R)
    for i = 2:k2-1
        facki_[i] = facki_[i-1]*ki_[i]
    end
    sign = -1
    s_ = Vector{elem_type(R)}(undef, k2)
    s_[1] = sign*facki_[k2-1]
    for i = 2:k2-1
        sign = -sign
        s_[i] = sign*facki_[i-1]*facki_[k2-i]
    end
    s_[k2] = facki_[k2-1]

    # ii. call algorithm 10_3
    points_ = elem_type(R)[ R(i) for i in 0:(k2-1) ]
    RPol = base_ring(M)
    x = gen(RPol)
    logk2 = floor(Int,log2(k2))
    Mij_ = Algorithm10_3(elem_type(RPol)[ x-points_[i] for i in 1:k2 ], logk2)

    # iii. call algorithm 10_9
    points_ = elem_type(R)[ R(i) for i in 0:(k2-1) ]
    Mk = zero(parent(M))
    for r = 1:n
        for c = 1:n
            Mk[r,c] = Algorithm10_9(points_,
                                    elem_type(R)[ s_[i]*MkEval_[i][r,c] for i in 1:k2 ], Mij_, RPol)
        end
    end

    return Mk
end

function MatrixEvaluation(Mk, k2, beta_)
    #    Implements [1, Corollary 10]
    #
    #    INPUT:
    #
    #    - ``Mk`` - a matrix of polynomials
    #    - ``k2`` - an integer >0
    #    - ``beta_``
    #
    #    OUTPUT:
    #    The sequence of matrixes M_k(beta) with beta in beta_
    #
    #    NOTE:
    #    If C = \ceil(beta_/k2) and M an nxn matrix, then the
    #    complexity is O( C(MM(n)k + n^2M(k)\log(k)) )
    #    (together w/ MatrixEvaluationPre)

    R = base_ring(base_ring(Mk))
    n = nrows(Mk)
    RPol = base_ring(Mk)
    x = gen(RPol)
    logk2 = floor(Int,log2(k2))

    # 3: Evaluate M_k(X) at the points specified by beta_
    res_ = dense_matrix_type(R)[zero_matrix(R, n, n) for i in 1:length(beta_)]
    p = 1
    cache = Dict{Int, Tuple{Vector{elem_type(RPol)}, Vector{elem_type(RPol)}}}()
    cache_10_3 = Dict{Int, Vector{Vector{elem_type(RPol)}}}()
    while (p+k2-1 <= length(beta_))
        moduli_ = elem_type(RPol)[ x-beta_[i] for i in p:(p+k2-1) ]
        if haskey(cache_10_3, logk2)
          _cache = cache_10_3[logk2]
        else
          _cache = Vector{elem_type(RPol)}[ elem_type(RPol)[ zero(RPol) for i in 1:2^(j - 1) ] for j in (logk2 + 1):-1:1]
          cache_10_3[logk2] = _cache
        end
        Mij_ = Algorithm10_3_cached(moduli_, logk2, _cache)
        if haskey(cache, length(moduli_))
          cache1, cache2 = cache[length(moduli_)]
        else
          cache1 = elem_type(RPol)[ zero(RPol) for i in 1:length(moduli_) ]
          cache2 = elem_type(RPol)[ zero(RPol) for i in 1:length(moduli_) ]
          cache[length(moduli_)] = (cache1, cache2)
        end
        for r = 1:n
            for c = 1:n
                values_ = Algorithm10_5_cached(Mk[r,c], moduli_, Mij_, cache1, cache2)
                #@assert values_ == Algorithm10_5(Mk[r,c], moduli_, Mij_)
                for i = 1:k2
                  res_[p+i-1][r,c] = deepcopy(values_[i])
                end
            end
        end
        p = p+k2
    end
    if (p <= length(beta_))
        b = length(beta_)
        for i = (b+1):(p+k2-1)
            push!(beta_, 0)
        end
        moduli_ = elem_type(RPol)[ x-beta_[i] for i in p:(p+k2-1) ]
        
        if haskey(cache_10_3, logk2)
          _cache = cache_10_3[logk2]
        else
          _cache = Vector{elem_type(RPol)}[ elem_type(RPol)[ zero(RPol) for i in 1:2^(j - 1) ] for j in (logk2 + 1):-1:1]
          cache_10_3[logk2] = _cache
        end

        Mij_ = Algorithm10_3_cached(moduli_, logk2, _cache)

        if haskey(cache, length(moduli_))
          cache1, cache2 = cache[length(moduli_)]
        else
          cache1 = elem_type(RPol)[ zero(RPol) for i in 1:length(moduli_) ]
          cache2 = elem_type(RPol)[ zero(RPol) for i in 1:length(moduli_) ]
          cache[length(moduli_)] = (cache1, cache2)
        end

        for r = 1:n
            for c = 1:n
                values_ = Algorithm10_5_cached(Mk[r,c], moduli_, Mij_, cache1, cache2)
                #@assert values_ == Algorithm10_5(Mk[r,c], moduli_, Mij_)
                for i = 1:b-p+1
                  res_[p+i-1][r,c] = deepcopy(values_[i])
                end
            end
        end
    end

    return res_
end

function UpperCaseDD(alpha, beta, k)

    #
    #    Implements [1, Lemma 6]
    #
    #    INPUT:
    #
    #    -  ``alpha`` - an element of R
    #    -  ``beta`` - an element of R
    #    -  ``k`` - a finite cardinal
    #
    #    OUTPUT:
    #
    #    Returns the element D(alpha, beta, k) of R

    k_ = Int[ k ]
    while k_[end] > 1
        push!(k_, floor(k_[end]/2))
    end

    R = parent(alpha)
    res = one(R)
    for i = 1:(length(k_)-1)
        res = res*LowerCaseDD(beta*(k_[i+1]+1),beta,k_[i+1])* LowerCaseDD(alpha*k_[i+1],beta,k_[i+1])
    end
    return res
end

function LinearRecurrence(M, L_, R_)
    s = Int(floor(log(4,R_[end])))
    R = base_ring(base_ring(M))
    return LinearRecurrence(M, L_, R_, inv(R(UpperCaseDD(1,2^s,2^s))), s)
end

function LinearRecurrence(M, L_,R_, DDi, s)
    #    Implements [2, Theorem 10]
    #
    #    INPUT:
    #
    #    -  ``M`` - a matrix of linear polynomials over R
    #    -  ``L_`` - a sequence of finite cardinals
    #    -  ``R_`` - a sequence of finite cardinals
    #    -  ``DDi`` - an (invertible) element of R
    #    -  ``s`` - a finite cardinal
    #
    #    OUTPUT:
    #
    #    A sequence of matrices of polynomials over R such that the i-th entry is
    #    given by M(L_[i],R_[i]) = M(L_[i]+1)M(L_[i]+2)...M(R_[i])
    #
    #    NOTE:
    #
    #    It is assumed that 0 <= L_[i] < R_[i] <= L_[i+1] < R_[i+1] <= T,
    #    where r = #L_ = #R_ and 0 < r < T^(1/2), s = Floor(log_4(R_[r]))
    #    and DDi = DD(1,2^s,2^s)^(-1)
    #
    #    The complexity is O( MM(n)T^(0.5) + n^2*M(T^(0.5)) ) (where M is an
    #    nxn matrix)

    # Check user input
    length(R_) != length(L_) && error("Number of interval boundaries do not match.")

    RPol = base_ring(M)
    x = gen(RPol)
    R = base_ring(RPol)

    # PRELIMINARIES
    r = length(L_)	# r is allowed to change (e.g. in the final step)
    r0 = length(L_)	# r0 must never change
    k = 2^s	# current interval length
    println(">> seeti")
    factorsDD_ = UpperCaseDD_(R(1),R(k),k)   # r_[1] = dd(alpha,beta,d)
    # with alpha = k*(d+1), beta=k, d=k/2
    # DDi_[1] = dd(alpha, beta, d)^(-1)
    DDi_ = RetrieveInverses(DDi,factorsDD_)
    println(">> invs1")
    if (k > 1)
        d = k >> 1
        # the 1st entry is 2^s
        factorsdd_ = LowerCaseDD_(R(k*(d+1)),R(k),d)

        # the entries d+1:3d+1 are 2^s,2*2^s, ...,(2^s+1)*2^s
        ki_ = RetrieveInverses(DDi_[1], factorsdd_)
        ki_ = elem_type(R)[ factorsdd_[1]*ki_[i] for i in (d+1):(length(factorsdd_)-1)]
        # ki_ = inverses of 1,2,...,k=2^s (in that order)
    else
        # ugly hack, in case k=1 we could do things much simpler
        ki_ = elem_type(R)[ R(1) ]
    end
    println(">> asdas")

    n = nrows(M)
    # res_[i] will contain M(L_i,R_i) at the end
    res_ = dense_matrix_type(R)[ identity_matrix(R, n) for j in 1:r0 ] # TODO this should probably be in Rmat????

    # STEP 0
    # get interval indices
    # l_[j] stores "index" of first interval
    # that fits into [L_[j]+1,...,R_[j]]
    l_ = Int[0 for j in 1:r0]
    for j = 1:r0
        qu,re = divrem(L_[j], k)
        if (re > 0)
            l_[j] = qu+2
        else
            l_[j] = qu+1
        end
    end
    # r_[j] stores "index" of last such interval
    r_ = Int[ div(R_[j],k) for j in 1:r0 ]

    # evaluate
    k_, logk, ddi1__, partiali1__, delta1__, s1_,
    ddi2__, partiali2__, delta2__, s2_ =
    MatrixAPEvaluationPre(k, R(1), R(k), DDi_, RPol)
    println(">> pre")
    # (the following currently accounts for 50% of all
    #  computation time in LinearRecurrences)
    M_ = MatrixAPEvaluation(M, k_, logk,
                            ddi1__, partiali1__, delta1__, s1_,
                            ddi2__, partiali2__, delta2__, s2_, 1, k, DDi_)
    println(">>eval")
    pop!(M_)
    append!(M_,
              MatrixAPEvaluation(Evaluate(M,x+k^2), k_, logk,
                                 ddi1__, partiali1__, delta1__, s1_,
                                 ddi2__, partiali2__, delta2__, s2_, 1, k, DDi_))
    pop!(M_)
    append!(M_,
              MatrixAPEvaluation(Evaluate(M,x+2*k^2), k_, logk,
                                 ddi1__, partiali1__, delta1__, s1_,
                                 ddi2__, partiali2__, delta2__, s2_, 1, k, DDi_))
    pop!(M_)
    append!(M_,
              MatrixAPEvaluation(Evaluate(M,x+3*k^2), k_, logk,
                                 ddi1__, partiali1__, delta1__, s1_,
                                 ddi2__, partiali2__, delta2__, s2_, 1, k, DDi_))
    # M_[i] contains M((i-1)k,ik), i=1,...,4k,4k+1 i.e.
    # the matrix belonging:the i-th interval [(i-1)k+1,ik]
    println(">> val")

    # collect/glue
    for j = 1:r0
        for m = l_[j]:r_[j]
            res_[j] = mul!(res_[j], res_[j], M_[m])
        end
    end
    # the interval [LApprox_[j]+1,...,RApprox_[j]]
    LApprox_ = Int[ (l_[j]-1)*k for j in 1:r0 ]
    # is the current approximation (and subinterval) of
    # [L_[j]+1,...,R_[j]]
    RApprox_ = Int[ r_[j]*k for j in 1:r0 ]
    # the matrix res_[j] stores
    # the corresponding matrix M(LApprox_[j],RApprox_[j])
    #  RApprox_ (and LApprox_) replaces N_j^{()}
    #  in the originial BGS07 proof

    # STEPS i
    # it may happen that after step 0, LApprox_[j] > RApprox_[j]
    # in that case, we set these values to a new common value between
    # L_[j] and R_[j]
    for j = 1:r0
        if (LApprox_[j] > RApprox_[j])
            newValue = L_[j] + (R_[j] - L_[j] >> 1)
            LApprox_[j] = newValue
            RApprox_[j] = newValue
        end
    end
    println(res_)

    while (k > 2*r0)
        # at the start of each loop we have
        # LApprox_[j]-k < L_[j] \le LApprox[j]
        # RApprox_[j] \le R_[j] < RApprox[j]+k
        kOld = k
        k = isqrt(r*k)
        if (k^2 < r*kOld)
            k = k+1
        end
        # the new interval length k = Ceiling(Sqrt(r*k))
        # note: since k > 2r, we have
        # k < kOld \le 2^s (*)
        # k2 is smallest power of 2 strictly larger than k
        k2 = 2^(floor(Int,log2(k))+1)
        # we need k2:use the tree structure of MatrixEvaluation
        # note: because of (*), k2 \le 2^s
        # for each interval boundary we'll have rho intervals
        # to better approximate it
        rho = div(kOld,k)

        # l_[j] = number 1,2,3,... of leftmost interval
        l_ = Int[ div((LApprox_[j] - L_[j]),k) for j in 1:r0 ]
        # (counting from the right and
        # starting from current approximation invertal)
        # of current interval length k that fits into\
        # [L_[j]+1,...,R_[j]]
        # likewise, r_[j] contains the number the rightmost interval
        # if l_[j] (resp. r_[j]) is 0,
        # then LApprox_[j]-L_[j]<k (resp. R_[j]-RApprox_[j]<k)
        # and there is nothing:do in this iteration
        r_ = Int[ div((R_[j] - RApprox_[j]), k) for j in 1:r0 ]

        # prepare MatrixEvaluation
        MPre = MatrixEvaluationPre(M, k, ki_[1:k2-1])
        # the ki_ elements that we need here
        # are all defined because of (*)

        # evaluate left
        evalPoints_ = Int[  ]
        for j = 1:r0
            evalPoints_ = append!(evalPoints_, Int[ LApprox_[j] - m*k for
                                                   m in 1:rho ])
        end
        M_ = MatrixEvaluation(MPre, k2, evalPoints_)
        # collect/glue left
        for j = 1:r0
            for m = 1:(l_[j])
                res_[j] = M_[m+(j-1)*rho]*res_[j]
            end
        end

        # evaluate right
        evalPoints_ = Int[ ]
        for j = 1:r0
            evalPoints_ = append!(evalPoints_, Int[ RApprox_[j] + m*k for
                                             m in 0:(rho-1) ])
        end
        M_ = MatrixEvaluation(MPre, k2, evalPoints_)
        # this could actually be optimised
        # in the first call:MatrixEvaluation in this loop
        # also output M_k(X)
        # and in this call simply evaluate
        # this wouldn't change overall O-complexity, though

        # collect/glue right
        for j = 1:r0
            for m = 1:r_[j]
                res_[j] = res_[j]*M_[m+(j-1)*rho]
            end
        end

        LApprox_ = Int[ LApprox_[j] - l_[j]*k for j in 1:r0 ]
        RApprox_ = Int[ RApprox_[j] + r_[j]*k for j in 1:r0 ]
    end
    println(res_)

    # FINAL STEP (k \le 2r)
    while (r >= 1)
        # at the start of each loop we have
        # LApprox_[j]-2*r < L_[j] \le LApprox[j]
        # RApprox_[j] \le R_[j] < RApprox[j]+2*r
        if (k > r)
            # else we may safely skip the following
            # also: we wouldn't have enough inverses
            r2 = 2^(floor(Int,log2(r))+1)
            # same thing as before with k and k2
            MPre = MatrixEvaluationPre(M, r, ki_[1:r2-1])
            M_ = MatrixEvaluation(MPre, r2, Int[ LApprox_[j]-r for j in 1:r0 ])
            for j = 1:r0
                if (LApprox_[j]-L_[j] >= r)
                    res_[j] = mul!(res_[j], M_[j], res_[j])
                    LApprox_[j] = LApprox_[j]-r
                end
            end
            M_ = MatrixEvaluation(MPre, r2, RApprox_)
            for j = 1:r0
                if (R_[j]-RApprox_[j] >= r)
                    res_[j] = mul!(res_[j], res_[j], M_[j])
                    RApprox_[j] = RApprox_[j]+r
                end
            end
            # now:
            # LApprox_[j] - L_[j] < r and
            # R_[j] - RApprox_[j] < r
        end
        if (r > 1)
            r = (r >> 1) + ModByPowerOf2(r,1)
        else
            r = 0
        end
    end
    println(res_)

    return res_
end

function mod!(f::PolyElem{T}, h, g) where {T}
   if length(g) == 0
      throw(DivideError())
   end

   set!(f, h)

   if length(h) >= length(g)
      b = lead(g)
      if !isone(b)
        g = inv(b)*g
      end
      c = base_ring(f)()
      while length(f) >= length(g)
         l = -lead(f)
         for i = 1:length(g) - 1
            c = mul!(c, coeff(g, i - 1), l)
            u = coeff(f, i + length(f) - length(g) - 1)
            u = addeq!(u, c)
            f = setcoeff!(f, i + length(f) - length(g) - 1, u)
         end
         set_length!(f, normalise(f, length(f) - 1))
      end
   end
   return f
end

function set!(a::PolyElem{T}, b::PolyElem{T}) where {T}
   fit!(a, length(b))
   for i = 1:length(b)
       a.coeffs[i] = set!(a.coeffs[i], b.coeffs[i])
   end
   set_length!(a, length(b))
   return a
end

function set!(z::padic, a::padic)
  ccall((:padic_set, :libflint), Nothing,
        (Ref{padic}, Ref{padic}, Ref{FlintPadicField}), z, a, parent(a))
  z.N = a.N
  return z
end

function set!(z::qadic, a::qadic)
   ccall((:qadic_set, :libflint), Nothing,
         (Ref{qadic}, Ref{qadic}, Ref{FlintQadicField}), z, a, parent(a))
   z.N = a.N
   return z
end

function set!(z::fmpq_abs_series, a::fmpq_abs_series)
   ccall((:fmpq_poly_init, :libflint), Nothing, (Ref{fmpq_abs_series},), z)
   ccall((:fmpq_poly_set, :libflint), Nothing,
        (Ref{fmpq_abs_series}, Ref{fmpq_abs_series}), z, a)
   z.prec = a.prec
   z.parent = parent(a)
   return z
end

