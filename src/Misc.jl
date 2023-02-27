
import AbstractAlgebra.Generic.LaurentSeriesElem
import AbstractAlgebra.Generic.LaurentSeriesFieldElem
import AbstractAlgebra.Generic.PolyElem


function Nemo.order(x::FinFieldElem)
    i = 1
    pow = x
    while pow != 1
        pow *= x
        i += 1
    end
    return i
end

#function gen(R::Nemo.GaloisField)
#    a = rand(R)
#    targ = modulus(R) - 1
#    while order(a) != targ
#        a = rand(R)
#    end
#    return a
#end
#Base.log(a::fmpz) = log(BigInt(a))

function (f::LaurentSeriesFieldElem)(x::LaurentSeriesFieldElem)
    @assert parent(f) == parent(x)
    K = parent(f)
    ret = zero(K)
    for i in 0:(precision(f) - 1)-valuation(f)
        ret += coeff(f, i+valuation(f))*(x^(i+valuation(f)))
    end
    return ret + O(gen(parent(f))^precision(f))
end



function Generic.mod(f::Generic.PolyElem{T}, g::Generic.PolyElem{T}) where {T <: Generic.SeriesElem}
   check_parent(f, g)
   if length(g) == 0
      throw(DivideError())
   end
   if length(f) >= length(g)
      f = deepcopy(f)
      b = leading_coefficient(g)
      g = inv(b)*g
      c = base_ring(f)()
      while length(f) >= length(g)
         l = -leading_coefficient(f)
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

