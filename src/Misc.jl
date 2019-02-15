

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
