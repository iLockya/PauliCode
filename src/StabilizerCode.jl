import Oscar: nvars, level, dim

abstract type AbstractLattice end
struct PBCLattice{n} <: AbstractLattice
    R∞ :: Union{MPolyRing, MPolyQuoRing}
    RL :: Union{MPolyRing, MPolyQuoRing}
    vars
    bdry
end

nvars(L::PBCLattice) = length(L.vars)
dim(L::PBCLattice) = nvars(L)
level(L::PBCLattice{n}) where n = n 
function Base.show(io::IO, L::PBCLattice)
    n = nvars(L)
    print(io, "Variables: ")
    for (i,v) in enumerate(L.vars)
        i != length(L.vars) ? print(io, "$v, ") : print(io, "$v.\n")
    end
    print(io, "Boundary conditions: ")
    if iszero(L.bdry)
        for (i,poly) in enumerate(L.R∞.I.gens)
            i != length(L.R∞.I.gens) ? print(io, "$poly, ") : print(io, "$poly.\n")
        end
    else
        for bi in eachcol(L.bdry)
            if ~iszero(bi)
                s = prod(L.R∞[j]^bi[j] for j in 1:n)
                print(io, "$s-1  ")
            end
        end
        print(io, ".\n")
    end
end

function oneDLattice(bdry = []; varlist::Vector{T} = [:y], n::Int = 2, kw...) where T<:Union{Char, Symbol, String}
    varlist = Symbol.(varlist)
    varlist_inv = cat(varlist, [Symbol(s,"̄") for s in varlist], dims=1)
end



"""
    SquareLattice(bdry = []; varlist::Vector{T} = [:x, :y], n::Int = 2, kw...) where T<:Union{Char, Symbol, String}
    
Construct a 2D square lattice with boundary conditions defined by `bdry`. Variables are `x` and `y` by default.

# Arguments
- `bdry::Vector{String}`: A vector of strings representing the boundary conditions.
- `varlist::Vector{T}`: A vector of symbols representing the variables. Default is `[:x, :y]`.
- `p::Int`: The prime number for the finite field. Default is `2`.

# Examples
```jldoctest
julia> L, (x,y) = SquareLattice([12,12])
julia> L.RL
Quotient
  of multivariate polynomial ring in 4 variables x, y, x̄, ȳ
    over prime field of characteristic 2
  by ideal (x*x̄ + 1, y*ȳ + 1, x^12 + 1, y^12 + 1)
```
"""
function SquareLattice(bdry = []; varlist::Vector{T} = [:x, :y], n::Int = 2, kw...) where T<:Union{Char, Symbol, String}
    varlist = Symbol.(varlist)
    varlist_inv = cat(varlist, [Symbol(s,"̄") for s in varlist], dims=1)

    if is_prime(n)
        R, xy = polynomial_ring(GF(n), varlist_inv; kw...)
    else
        R, xy = polynomial_ring(residue_ring(ZZ,n)[1], varlist_inv; kw...)
    end

    I∞ = ideal(R, [xy[1]*xy[3]-1, xy[2]*xy[4]-1])
    R∞, p0 = quo(R, I∞)
    eval(:($(Symbol(varlist[1])) = $(p0(xy[1]))))
    eval(:($(Symbol(varlist[2])) = $(p0(xy[2]))))

    if isempty(bdry)
        bdry_matrix = [0 0; 0 0]
    elseif bdry isa Vector{Int}
        @assert length(bdry) == 2
        bdry_matrix = [bdry[1] 0; 0 bdry[2]]
    elseif bdry isa Matrix
        @assert size(bdry) == (2, 2)
        bdry_matrix = bdry
    else
        error("Boundary conditions must be a vector of integers or a 2x2 matrix.")
    end
    Ib = ideal(R∞, [R∞[1]^bdry_matrix[1,1]*R∞[2]^bdry_matrix[1,2] - 1,
                    R∞[1]^bdry_matrix[2,1]*R∞[2]^bdry_matrix[2,2] - 1])
    RL, _ = quo(R∞, Ib)
    xy = RL.(xy)
    return PBCLattice{n}(R∞, RL, xy[1:2], bdry_matrix), xy[1:2]
end

"""
    CubicLattice(bdry = []; varlist::Vector{T} = [:x, :y, :z], n::Int = 2, kw...) where T<:Union{Char, Symbol, String}
    
Construct a 3D square lattice with boundary conditions defined by `bdry`. Variables are `x`, `y` and `z` by default.

# Arguments
- `bdry::Vector{String}`: A vector of strings representing the boundary conditions.
- `varlist::Vector{T}`: A vector of symbols representing the variables. Default is `[:x, :y]`.
- `p::Int`: The prime number for the finite field. Default is `2`.

# Examples
```jldoctest
julia> L, (x,y,z) = CubicLattice([12,12,12])
julia> L.RL
Quotient
  of multivariate polynomial ring in 6 variables x, y, z, x̄, ȳ, z̄
    over prime field of characteristic 2
  by ideal (x*x̄ + 1, y*ȳ + 1, z*z̄ + 1, x^12 + 1, y^12 + 1, z^12 + 1)
```
"""
function CubicLattice(bdry = []; varlist::Vector{T} = [:x, :y, :z], n::Int = 2, kw...) where T<:Union{Char, Symbol, String}
    varlist = Symbol.(varlist)
    varlist_inv = cat(varlist, [Symbol(s,"̄") for s in varlist], dims=1)

    if is_prime(n)
        R, xyz = polynomial_ring(GF(n), varlist_inv; kw...)
    else
        R, xyz = polynomial_ring(residue_ring(ZZ,n)[1], varlist_inv; kw...)
    end

    I∞ = ideal(R, [xyz[1]*xyz[4]-1, xyz[2]*xyz[5]-1, xyz[3]*xyz[6]-1])
    R∞, p0 = quo(R, I∞)
    eval(:($(Symbol(varlist[1])) = $(p0(xyz[1]))))
    eval(:($(Symbol(varlist[2])) = $(p0(xyz[2]))))
    eval(:($(Symbol(varlist[3])) = $(p0(xyz[3]))))

    if isempty(bdry)
        bdry_matrix = zeros(Int, 3,3)
    elseif bdry isa Vector{Int}
        @assert length(bdry) == 3
        bdry_matrix = [bdry[1] 0 0; 0 bdry[2] 0; 0 0 bdry[3]]
    elseif bdry isa Matrix
        @assert size(bdry) == (3, 3)
        bdry_matrix = bdry
    else
        error("Boundary conditions must be a vector of integers or a 3x3 matrix.")
    end

    Ib = ideal(R∞, [R∞[1]^bdry_matrix[1,1]*R∞[2]^bdry_matrix[1,2]*R∞[3]^bdry_matrix[1,3] - 1,
                    R∞[1]^bdry_matrix[2,1]*R∞[2]^bdry_matrix[2,2]*R∞[3]^bdry_matrix[2,3] - 1,
                    R∞[1]^bdry_matrix[3,1]*R∞[2]^bdry_matrix[3,2]*R∞[3]^bdry_matrix[3,3] - 1])
    RL, _ = quo(R∞, Ib)
    xyz = RL.(xyz)
    return PBCLattice{n}(R∞, RL, xyz[1:3], bdry_matrix), xyz[1:3]
end

struct ClassicalCode{n} <: AbstractStabilizerCode{n}
    lattice::AbstractLattice
    X::FreeModuleHom
    function ClassicalCode{n}(lattice::AbstractLattice, X::FreeModuleHom) where n
        matrix(X)
        image(X)
        new{n}(lattice, X)
    end
end

struct CSSCode{n} <: AbstractStabilizerCode{n}
    lattice::AbstractLattice
    X::FreeModuleHom
    Z::FreeModuleHom
    function CSSCode{n}(lattice::AbstractLattice, X::FreeModuleHom, Z::FreeModuleHom) where n
        matrix(X), matrix(Z)
        image(X), image(Z)
        new{n}(lattice, X, Z)
    end
end

function Base.show(io::IO, C::CSSCode)
    print(io, "Variables: ")
    for (i,v) in enumerate(C.lattice.vars)
        i != nvars(C) ? print(io, "$v, ") : print(io, "$v.\n")
    end
    println(io, "𝕏 checks:")
    display(transpose(matrix(C.X)))
    println(io, "ℤ checks:")
    display(transpose(matrix(C.Z)))
    print(io, "Boundary conditions: ")
    if iszero(C.lattice.bdry)
        for (i,poly) in enumerate(C.lattice.R∞.I.gens)
            i != length(C.lattice.R∞.I.gens) ? print(io, "$poly, ") : print(io, "$poly ].\n")
        end
    else
        for bi in eachcol(C.lattice.bdry)
            if ~iszero(bi)
                s = prod(C.lattice.R∞[j]^bi[j] for j in 1:nvars(C))
                print(io, "$s-1 ")
            end
        end
        print(io, ".\n")
    end
end

nvars(C::AbstractStabilizerCode) = nvars(C.lattice)
level(C::AbstractStabilizerCode{n}) where n = n

"""
    ClassicalCode(X_checks::VecOrMat{String}, bdry::Vector{String}=[]; p::Int=2, kw...)

Construct a classical stabilizer Hamiltonian (X check only).

# Arguments
- `X_checks::VecOrMat{String}`: `n×q` matrix of X checks, where `n` is the number of checks and `q` is the number of qubits per unit cell.
- `bdry::Vector{String}`: Boundary condition.
- `p::Int`: qubit/qudit level.

# Examples
```jldoctest
julia> C, (x,y) = ClassicalCode(["1+x" "0";"0" "1+y"], ["x^10-1", "y^10-1"]); # 2D Ising model.
julia> C.Z.matrix
[x + 1       0]
[    0   y + 1]
```
"""
function ClassicalCode(X_checks::VecOrMat{String}, bdry::VecOrMat=[]; p::Int=2, kw...)
    if ndims(X_checks) == 1
        X_checks = reshape(X_checks, (1, :))
    end
    varlist = extract_variables([X_checks..., bdry...])
    if length(varlist) == 2
        lattice, xyz = SquareLattice(bdry, varlist=varlist, p=p; kw...)
        eval(:($(varlist[1]) = $xyz[1]))
        eval(:($(varlist[2]) = $xyz[2]))
    elseif length(varlist) == 3
        lattice, xyz = CubicLattice(bdry, varlist=varlist, p=p; kw...)
        eval(:($(varlist[1]) = $xyz[1]))
        eval(:($(varlist[2]) = $xyz[2]))
        eval(:($(varlist[3]) = $xyz[3]))
    else
        error("Only 2D and 3D lattices are supported.")
    end
    n, q = size(X_checks)
    F = free_module(lattice.R∞, n)
    G = free_module(lattice.R∞, q)
    V = typeof(G[1])[]
    for i in 1:n
        check = zero(G)
        for j in 1:q
            check += eval(Meta.parse(X_checks[i, j])) * G[j]
        end
        push!(V, check)
    end
    XStab = hom(F, G, V)
    return ClassicalCode{p}(lattice, XStab), xyz
end

"""
    CSSCode(X_checks::T, Z_checks::T, bdry::Vector{String}; n::Int=2, kw...) where T<:VecOrMat{String}
    CSSCode(X_checks::T, bdry = []; p::Int=2, kw...) where T<:VecOrMat{String}

Construct a quantum stabilizer Hamiltonian.

If `Z_checks` is not provided, it is assumed to be the hermitian conjugate of `X_checks`.

# Arguments
- `X_checks::VecOrMat{String}`: `t×q` matrix of X checks, where `t` is the number of checks and `q` is the number of qubits per unit cell.
- `Z_checks::VecOrMat{String}`: Same as `X_checks`.
- `bdry::Vector{String}`: Boundary condition.
- `p::Int`: qubit/qudit level.
# Examples
```jldoctest
julia> C, (x,y) = CSSCode(["1+x", "1+y"]) # Toric code.
julia> C.X.matrix
[x + 1   y + 1]
```
"""
function CSSCode(X_checks::VecOrMat{T}, Z_checks::VecOrMat{T}, bdry::VecOrMat{Int}; n::Int=2, kw...) where T<:String
    if ndims(X_checks) == 1
        X_checks = reshape(X_checks, (1, :))
    end
    if ndims(Z_checks) == 1
        Z_checks = reshape(Z_checks, (1, :))
    end
    varlist = extract_variables([X_checks..., Z_checks...])
    if length(varlist) == 2
        lattice, xyz = SquareLattice(bdry, varlist=varlist, n=n; kw...)
        eval(:($(varlist[1]) = $xyz[1]))
        eval(:($(varlist[2]) = $xyz[2]))
    elseif length(varlist) == 3
        lattice, xyz = CubicLattice(bdry, varlist=varlist, n=n; kw...)
        eval(:($(varlist[1]) = $xyz[1]))
        eval(:($(varlist[2]) = $xyz[2]))
        eval(:($(varlist[3]) = $xyz[3]))
    else
        error("Only 2D and 3D lattices are supported.")
    end
    t, q = size(X_checks)
    G = free_module(lattice.RL, q)
    F = free_module(lattice.RL, t)
    V = typeof(G[1])[]
    for i in 1:t
        check = zero(G)
        for j in 1:q
            check += eval(Meta.parse(X_checks[i, j])) * G[j]
            println(Meta.parse(X_checks[i, j]))
        end
        push!(V, check)
    end
    XStab = hom(F, G, V)
    t, q = size(Z_checks)
    F = free_module(lattice.RL, t)
    V = typeof(G[1])[]
    for i in 1:t
        check = zero(G)
        for j in 1:q
            check += eval(Meta.parse(Z_checks[i, j])) * G[j]
        end
        push!(V, check)
    end
    ZStab = hom(F, G, V)
    return CSSCode{n}(lattice, XStab, ZStab), xyz
end

CSSCode(X_checks::VecOrMat{String}, bdry::VecOrMat; n::Int = 2, kw...) = 
  CSSCode(X_checks, antipode.(reverse(X_checks)), bdry, n=n; kw...)
CSSCode(X_checks::VecOrMat{String}; n::Int = 2, kw...) = CSSCode(X_checks, Int[]; n=n, kw...)


#_x_excitation_map(φ::FreeModuleHom) = hom(codomain(φ), domain(φ), transpose_with_antipode(matrix(φ)))
#_z_excitation_map(φ::FreeModuleHom) = hom(codomain(φ), domain(φ), -transpose_with_antipode(matrix(φ)))
_z_excitation_map(C::CSSCode) = hom(codomain(C.X), domain(C.X), transpose_with_antipode(matrix(C.X)))
_x_excitation_map(C::CSSCode) = hom(codomain(C.Z), domain(C.Z), -transpose_with_antipode(matrix(C.Z)))
"""
    excitation_map(C::CSSCode)
"""
function excitation_map(C::CSSCode)
    _x_excitation_map(C), _z_excitation_map(C)
end


function symplectic_kernel(M)
    T = eltype(M)
    if T <: LaurentPolyRingElem
        return antipode.((kernel(M)))
    else
        return antipode.(base_ring(M).(kernel(lift.(M))))
    end
end