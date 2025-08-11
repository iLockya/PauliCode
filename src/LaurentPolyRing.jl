abstract type AbstractLauPolyRing <: Ring end
abstract type AbstractLauPolyRingElem <: RingElem end

struct LauPolyRing{T} <: AbstractLauPolyRing
    R :: MPolyQuoRing{T}
end
function Base.show(io::IO, R::LauPolyRing)
    println(io, "Univariate Laurent polynomial ring  in $(gens(R)[1]) over $(coefficient_ring(R)).")
end

struct LauMPolyRing{T} <: AbstractLauPolyRing
    R :: MPolyQuoRing{T}
end
function Base.show(io::IO, R::LauMPolyRing)
    println(io, "Multivariate Laurent polynomial ring in $(gens(R))")
    println(io, "  over $(coefficient_ring(R)).")
end

##### methods with MPolyQuoRing #####
Base.one(R::AbstractLauPolyRing) = one(R.R)
Base.zero(R::AbstractLauPolyRing) = zero(R.R)
Oscar.coefficient_ring(R::AbstractLauPolyRing) = coefficient_ring(R.R)
Oscar.base_ring(R::AbstractLauPolyRing) = coefficient_ring(R.R)
Oscar.nvars(R::AbstractLauPolyRing) = nvars(base_ring(R.R))÷2
Oscar.gens(R::AbstractLauPolyRing) = first(gens(base_ring(R.R)), nvars(R))
Oscar.gen(R::AbstractLauPolyRing, i::Int) = i <= nvars(R) ? gens(R)[i] : "$i out of range."
Oscar.free_module(R::AbstractLauPolyRing, n::Int) = free_module(R.R, n)

struct LauPolyRingElem{T} <: AbstractLauPolyRingElem
    x :: MPolyQuoRingElem{T}
    parent :: R :: LauPolyRing{T}
end

struct LauMPolyRingElem{T} <: AbstractLauPolyRingElem
    x :: MPolyQuoRingElem{T}
    parent :: R :: LauMPolyRing{T}
end
Base.show(io::IO, x::AbstractLauPolyRingElem) = print(io, x.x)

function laurent_poly_ring(R::Ring, s::Union{Char, AbstractString, Symbol}; kw...)
    v = string(s)
    R, vs = polynomial_ring(R, [v, v*"̄"]; kw...)
    A, _ = quo(R, ideal(prod(vs)-1))
    return LauPolyRing(A), LauPolyRingElem(A(vs[1]))
end

function laurent_poly_ring(R::Ring, ss::Vector{Symbol}; kw...)
    v = string.(ss)
    n = length(v)
    v_inv = [x * "̄" for x in v]
    R, vs = polynomial_ring(R, cat(v, v_inv, dims=1); kw...)
    I = ideal(R, [vs[i]*vs[n+i] - 1 for i in 1:n])
    A, _ = quo(R, I)
    return LauMPolyRing(A), LauMPolyRingElem.(A.(vs[1:n]))
end


