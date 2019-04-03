
import AbstractAlgebra.Generic.LaurentSeriesElem

#@doc Markdown.doc"""
#    root(a::fmpq, n::Int)
#> Given $x$ and in integer $n$, returns $\sqrt[n]{x}$ if it exists.
#"""
function Nemo.root(x::fmpq, n::Integer)
   nnum = root(numerator(x), n)
   (nnum^n != numerator(x)) && throw(DomainError("Argument $x has no $n th root"))
   nden = root(denominator(x), n)
   (nden^n != denominator(x)) && throw(DomainError("Argument $x has no $n th root"))
   return nnum//nden
end


function (R::Nemo.GaloisField)(x::fmpq)
    return R(numerator(x))//R(denominator(x))
end


function Nemo.add!(a::padic, b::padic, c::padic)
   return b + c
end

function (f::LaurentSeriesFieldElem)(x::LaurentSeriesFieldElem)
    @assert parent(f) == parent(x)
    K = parent(f)
    ret = zero(K)
    for i in 0:(precision(f) - 1)-valuation(f)
        ret += coeff(f, i+valuation(f))*(x^(i+valuation(f)))
    end
    return ret + O(gen(parent(f))^precision(f))
end


function Generic.integral(x::LaurentSeriesElem)
   z = deepcopy(x)
   set_prec!(z, precision(x) + 1)
   set_val!(z, valuation(x) + 1)
   len = pol_length(x)
   for i = 1:len
       z = setcoeff!(z, i - 1, polcoeff(x, i - 1)//((i-1)*Generic.scale(z)+valuation(z)))
   end
   set_length!(z, normalise(z, len))
   renormalize!(z)
   z = rescale!(z)
   return z
end

function Generic.derivative(x::LaurentSeriesElem{T}) where {T <: RingElement}
   z = deepcopy(x)
   set_prec!(z, precision(x) - 1)
   set_val!(z, valuation(x) - 1)
   len = pol_length(x)
   for i = 1:len
       z = setcoeff!(z, i - 1, ((i-1)*Generic.scale(z)+valuation(z) + 1)*polcoeff(x, i - 1))
   end
   set_length!(z, normalise(z, len))
   renormalize!(z)
   z = rescale!(z)
   return z
end
