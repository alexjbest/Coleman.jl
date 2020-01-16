#  Copyright (C) 2019 Alex J. Best <alex.j.best@gmail.com>
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
# This module implements the accumulating remainder tree of Costa, Gerbicz, Harvey
# [1] over arbitrary rings.
#
# [1] Costa, Edgar, Robert Gerbicz, and David Harvey. “A Search for Wilson
#     Primes.” Mathematics of Computation 83, no. 290 (2014): 3071–91.
#     https://doi.org/10.1090/S0025-5718-2014-02800-7.
#
#*****************************************************************************

include("LinearRecurrence.jl")
include("Misc.jl")

import AbstractAlgebra.RelSeriesElem
import AbstractAlgebra.Ring
import AbstractAlgebra.Generic
import AbstractAlgebra.Generic.LaurentSeriesElem
import AbstractAlgebra.Generic.LaurentSeriesFieldElem
using Nemo, Hecke

export ART


@doc doc"""
    get_primes(M, N)
> Return the list of all primes M < p <= N.
"""
function get_primes(M, N)
    # TODO ideally this should seive!
    return Hecke.PrimesSet(UInt(M), UInt(N))
end

@doc doc"""
    ART(n, lambda, B, M)
> Run the accumulating remainder tree for a set of $B$ total $n\times n$
> matrices $M(k)$ to compute the sequence of products
> $M_0 M_1 \cdots M_{(p-1)/2} \pmod {p^\lambda}$
> for all $p < 2B$.
> $M$ should be a function that returns the $k$th matrix in the list, but could just be a list accessor like `x -> getindex(M, x)`.
"""
function ART(n, lambda, B, M)
    primes = get_primes(2, 2*B + 1)
    # trees of depth l
    l = clog2(B)
    P = [[nothing for j in 1:2^i] for i in 1:l]
    for i in (l-1):-1:1
        for j in 1:2^i
            P[i][j] = P[i+1][2*j] * P[i+1][2*j+1]


        end
    end

    return
end
