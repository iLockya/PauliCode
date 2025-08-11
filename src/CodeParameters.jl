"""
    physical_qubits(L::AbstractLattice)
    physical_qubits(C::AbstractStabilizerHamiltonian)

Return the number of physical qubits in the code.

# Examples
```jldoctest
julia> C, (x,y) = CSSCode(["x^3+y^2+y","y^3+x^2+x"],["x^12-1","y^12-1"])
julia> physical_qubits(C)
144
```
"""
function physical_qubits(L::AbstractLattice)
    is_finite_dimensional_vector_space(L.RL) ? dim(L) * vector_space_dimension(L.RL) : Inf
end
physical_qubits(C::AbstractStabilizerCode) = physical_qubits(C.lattice)

function _lift_cokernel(φ::FreeModuleHom)
    img = sub_object(codomain(φ), φ.imgs_of_gens)
    lift_img = Oscar._poly_module(img.sum)
    quo_object(lift_img.F, gens(lift_img))
end

"""
    vdim_cokernel(φ::FreeModuleHom)

Return `vdim(cokerφ)`, force `domain(φ)` to be `Rb`. If `Rb` is not provided, return `vdim(cokerφ)`.

# Examples
```jldoctest
julia> Q, (x,y) = CSSCode(["x^3+y^2+y"; "y^3+x^2+x"],["x^10-1","y^10-1"]);
julia> vdim_cokernel_over_qring(Q.X)
0
julia> vdim_cokernel_over_qring(Q.X, Q.lattice.RL)
8
```
"""
function vdim_cokernel(φ::FreeModuleHom)
    M = _lift_cokernel(φ)
    vdim = vector_space_dimension(M)
    return vdim < 0 ? Inf : vdim
end

function vbase_cokernel(φ::FreeModuleHom)
    M = _lift_cokernel(φ)
    vector_space_basis(M)
end

    
    
"""
    topological_index(C::CSSCode{p}) where p

Return the topological index `Q` of the code. Equivalent to `GSD` at thermodynamic limit.

Reference: https://arxiv.org/abs/2503.04699
# Examples
```jldoctest
julia> C, (x,y) = CSSCode(["x^3+y^2+y","y^3+x^2+x"],[12,12])
julia> topological_index(C)
8
```
"""
function topological_index(C::CSSCode{p}) where p
    R = C.lattice.R∞
    (nx, q) = size(matrix(C.X))
    ϵX = hom(free_module(R, q), free_module(R, nx), 
           R.(transpose_with_antipode(matrix(C.X))))
    vdim_cokernel(ϵX)
end


"""
    GSD(C::CSSCode{p}) where p

Return the ground state degeneracy of the code Hamiltonian.
"""
function GSD(C::CSSCode{p}) where p
    ϵX, ϵZ = excitation_map(C)
    dx = _homological_dimension(C.X, ϵX)
    dz = _homological_dimension(C.Z, ϵZ)
    (dx + dz) ÷ 2
end

function fractal_charge(C::CSSCode{p}) where p
    ϵX, ϵZ = excitation_map(C)
    ϵX, ϵZ = change_base_ring(C.lattice.R∞, ϵX)[1], change_base_ring(C.lattice.R∞, ϵZ)[1]
    vdim_cokernel(ϵX) + vdim_cokernel(ϵZ)
end

function _homological_dimension(∂::T, ϵ::T) where T<:FreeModuleHom
    K = quo_object(kernel(ϵ)[1], image(∂)[1]) |> present_as_cokernel
    L = Oscar._poly_module(K.quo)
    Q = quo_object(L.F, gens(L))
    return vector_space_dimension(Q)
end
"""
    logical_qubits(C::CSSCode{p}) where p
    logical_qubits(C::ClassicalCode{p}) where p

Return the number of logical qubits in the code.
# Examples
```jldoctest
julia> C, (x,y) = CSSCode(["x^3+y^2+y","y^3+x^2+x"],[12,12])
julia> logical_qubits(C)
16
```
"""
function logical_qubits(C::CSSCode{p}) where p
    GSD(C)
end

function logical_qubits(C::ClassicalCode{p}) where p
    vdim_cokernel(C.X) # right?
end



struct MobilitySubLattice
    periods :: Dict
end
function Base.show(io::IO, M::MobilitySubLattice)
    for (var, period) in M.periods
        println(io, "  $var => $period")
    end
end

"""
    MobilitySubLattice(C::CSSCode)

Mobility sub lattice of a CSS code `C` is defined as the lattice of periods of the generators of the code.

# Example
```jldoctest
julia> C, (x,y) = CSSCode(["1+x+x^(-1)*y^(-3)","1+y+y^(-1)*x^3"],["x^12-1","y^12-1"]);
julia> MobilityLattice(C)
MobilitySubLattice(Integer[762, 762])
```
"""
function MobilitySubLattice(C::CSSCode)
    I = ideal(lift.([C.X.matrix...]))+C.lattice.R∞.I
    R = base_ring(I)
    n = nvars(C)
    generators = gens(R)
    periods = Dict{Symbol, Integer}()
    for i in 1:n
        ordering = lex([generators[1:i-1]; generators[i+1:end]; generators[i]])
        gb = standard_basis(I, ordering=ordering, complete_reduction=true)
        periods[Symbol(generators[i])] = polynomial_period(gb[1])
    end
    return MobilitySubLattice(periods)
end