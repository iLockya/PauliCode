
struct BoundaryStabilizer{T}
    origin :: T
    bulk :: T
    truncated :: T
end

abstract type AbstractBoundary end
abstract type AbstractCSSBoundary <: AbstractBoundary end
"""
    SmoothCSSBoundary{n} <: AbstractOpenBoundary

Truncate the x⁻ⁿ terms.
"""
struct SmoothCSSBoundary{n} <: AbstractCSSBoundary
    XB :: BoundaryStabilizer
    ZB :: BoundaryStabilizer
end

struct RoughCSSBoundary{n} <: AbstractCSSBoundary
    XB :: BoundaryStabilizer
    ZB :: BoundaryStabilizer
end


function _position_to_modulehom(R, P, d)
    y = gen(R, 1)
    F = free_module(R, length(P))
    G = free_module(R, d)
    φ = zeros(R, rank(F), rank(G))
    for (i, Pi) in enumerate(P)
        for (a, b) in Pi
            φ[i, a] += y^b
        end
    end
    hom(F, G, matrix(φ))
end

function smooth_boundary(C::CSSCode{n}) where n
    PX = _position_info_smooth(C.X)
    PZ = _position_info_smooth(C.Z)
    d = maximum([p[1] for p in cat(PX...,PZ...,dims=1)])
    if is_prime(n)
        # R, y = laurent_polynomial_ring(GF(n), :y)
        R, (y,ȳ) = polynomial_ring(GF(n), [:y,:ȳ])
        A, _  = quo(R, ideal(y*ȳ-1))
        y, ȳ = gens(A)
    else
        R, (y,) = polynomial_ring(residue_ring(ZZ,n)[1], [:y])
        R, y = laurent_polynomial_ring(residue_ring(ZZ,n)[1], :y)
    end
    # X boundary stabilizers
    PB, PT = _bulk_stabilizers(PX, d), _truncated_stabilizers(PX, d)
    XB = BoundaryStabilizer(_position_to_modulehom(A, PX, d),
                            _position_to_modulehom(A, PB, d),
                            _position_to_modulehom(A, PT, d))
    # Z boundary stabilizers
    PB, PT = _bulk_stabilizers(PZ, d), _truncated_stabilizers(PZ, d)
    ZB = BoundaryStabilizer(_position_to_modulehom(A, PZ, d),
                            _position_to_modulehom(A, PB, d),
                            _position_to_modulehom(A, PT, d))
    SmoothCSSBoundary{n}(XB, ZB)
end


function _bulk_stabilizers(P, d::Int)
    Sb = []
    for Pi in P
        Sbi = []
        push!(Sbi, Pi)
        while true
            S_k = Vector{Int}[]
                for (a, b) in Pi
                    a + 2 ≤ d && push!(S_k, [a + 2, b])
                end
            Pi = S_k
            isempty(Pi) ? break : push!(Sbi, Pi)
        end
        Sb = [Sb; Sbi]
    end
    Sb
end

function _truncated_stabilizers(P, d::Int)
    St = _bulk_stabilizers(P, d)
    for Pi in P
        Sti = []
        while true
            S_k = Vector{Int}[]
            for (a, b) in Pi
                a - 2 ≥ 1 && push!(S_k, [a - 2, b])
            end
            Pi = S_k
            isempty(Pi) ? break : push!(Sti, Pi)
    end
    St = [St; Sti]
    end
    St
end
"""
    _smooth_pos_map(V::Vector{Int}, j::Int) :: Vector{Int}

If `j==1`, `(a,b) ↦ (2a+2, b)`, if `j==2`, `(a,b) ↦ (2a+1, b)`.
"""
function _smooth_pos_map(V::Vector{Int}, j::Int)
    if j == 1
        return Vector{Int}([2*V[1]+2, V[2]])
    elseif j == 2
        return Vector{Int}([2*V[1]+1, V[2]])
    end
end

function _position_info_smooth(φ::FreeModuleHom)
    P = []
    t, q = size(matrix(φ))
    @assert q == 2 # two dimensions 
    for i in 1:t # enumerate checks
        M = φ.matrix[i,:]
        Pi = Vector{Int}[]
        max_inverse_pow = _max_inverse_pow(M)
        for (j,f) in enumerate(M) # horizontal/vertical checks
            for pows in exponents(f.f) # each checked qubits
                shift_pows = first(pows, 2) + max_inverse_pow - last(pows, 2)
                push!(Pi, _smooth_pos_map(shift_pows, mod1(j,2)))
            end
        end
        push!(P,Pi)
    end
    return P
end

function _max_inverse_pow(V)
    R = base_ring(V[1])
    q = nvars(R) ÷ 2
    max_inverse_pow = repeat(Int[0], q)
    for f in V
        for pows in exponents(f.f)
            max_inverse_pow = max.(max_inverse_pow, last(pows, q))
        end
    end
    return max_inverse_pow
end




function rough_boundary(C::CSSCode{n}) where n
end

function number_of_Z_excitations(B::AbstractCSSBoundary)
    ϵz = _z_excitation_map(B.ZB.truncated)
    vdim_cokernel(ϵz)
end
