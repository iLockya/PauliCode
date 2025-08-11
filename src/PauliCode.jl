module PauliCode


using Oscar 

abstract type AbstractStabilizerCode{p} end

include("StabilizerCode.jl")
include("CodeParameters.jl")
include("OpenBoundary.jl")
# include("LaurentPolyRing.jl")
include("Binomial/Binomial.jl")
include("utils.jl")


export SquareLattice, CubicLattice, ClassicalCode, CSSCode
export physical_qubits, logical_qubits, topological_index, GSD, vdim_cokernel,
        vbase_cokernel
export minimal_polynomial, polynomial_period, nvars, level, symplectic_kernel

export MobilitySubLattice

export Lattice, exponent_lattice_zero_dim

# Test functions
# export separable_decomposition_perfect # 
# export next_combination, has_next_combination,unipotent_lattice_p, primitive_element,
# is_primitive,exponent_lattice_finite_field

end
