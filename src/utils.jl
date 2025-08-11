"""
    irreducible_polynomial_period(f) :: Integer

Returns the period of the irreducible polynomial `f`. 
Period is defined as the smallest integer `n` such that `x^n = 1 mod f`.
# Example
```jldoctest
julia> S, x = polynomial_ring(GF(3), :x)
(Univariate polynomial ring in a over GF(3), x)
julia> polynomial_period(x^2+x+2)
8
```
"""
function irreducible_polynomial_period(f::Union{FqPolyRingElem, FqMPolyRingElem}) :: Integer
    @assert is_irreducible(f) "$f must be irreducible."
    @assert characteristic(f.parent) > 0 "$(f.parent) must be a finite field."
    if f isa FqMPolyRingElem
        @assert is_univariate(f) "f must be univariate."
        idx = argmax(degrees(f))
        a = f.parent[idx]
        m = degree(f, a)
    elseif f isa FqPolyRingElem
        a = gen(f.parent)
        m = degree(f)
    else
        error("Unsupported polynomial type.")
    end
    q = order(f.parent.base_ring)
    fac = factor(q^m-1).fac
    n = 1
    for (p, e) in fac
        for i in 1:e
            if mod(a^((q^m-1)÷p^i), f) != 1
                n *= p^(e+1-i)
                break
            end
        end
    end
    n
end

"""
    polynomial_period(f::Union{FqPolyRingElem, FqMPolyRingElem}) :: Integer

Returns the period of the polynomial `f`.
Period is defined as the smallest integer `n` such that `x^n = 1 mod f`.

# Example
```jldoctest
julia> R, x = polynomial_ring(GF(2), :x)
(Univariate polynomial ring in x over GF(2), x)
julia> polynomial_period(x^11+x^10+x^9+x^7+x^6+x^4+1)
762
```
"""
function polynomial_period(f::Union{FqPolyRingElem, FqMPolyRingElem}) :: Integer
    f isa FqMPolyRingElem && @assert is_univariate(f) "f must be univariate."
    p = characteristic(f.parent)
    fac = factor(f).fac
    power_of_p = 1
    n = 1
    for (fi, mi) in fac
        n *= irreducible_polynomial_period(fi)
        while power_of_p < mi
            power_of_p *= p
        end
    end
    n * power_of_p
end

function extract_variables(poly_array::Array{String})
    varlist = Set{String}()
    for poly in poly_array
        for m in eachmatch(r"[a-zA-Z]+[0-9]*", poly)
            push!(varlist, m.match)
        end
    end
    return Symbol.(collect(varlist))
end


"""
    antipode(f::MPolyQuoRingElem) :: MPolyQuoRingElem

Returns the antipode of the polynomial `f`. 
The antipode map is  `xᵃyᵇ ↦ x⁻ᵃy⁻ᵇ`.
"""
function antipode(f::MPolyQuoRingElem)  
    g = zero(f)
    for (mon, coeff) in zip(monomials(f.f), coefficients(f.f))
        g += coeff * f.P(mon)^(-1)
    end
    return g
end

function antipode(f::LaurentPolyRingElem)
    g = zero(f)
    g.mindeg = -(degree(f.poly) + f.mindeg)
    g.poly = reverse(f.poly)
    g
end
function antipode(f::String) :: String
    return replace(f, r"([a-zA-Z]+[0-9]*)" => s"(\1^-1)")
end

transpose_with_antipode(M) = transpose(antipode.(M))



