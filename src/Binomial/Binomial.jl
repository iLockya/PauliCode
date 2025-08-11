include("exponent_lattice_zero_dim.jl")
include("exponent_lattice_perffect_fields.jl")
include("exponent_lattice_extension_field.jl")
include("utils.jl")

export reduction_to_zero_dim, exponent_lattice, exponent_lattice_finite_field_max

function reduction_to_zero_dim(I, elems)
    R = base_ring(I)
    @assert I : ideal(elems) == I "`I` is not saturated wrt `elems`."
    dim(I) == 0 && return [I]
    # Compute independent sets by Sigular functions.
    Rs = singular_poly_ring(R)
    Is = Singular.Ideal(Rs, Rs.(gens(I))) |> Singular.std
    indep_indices = [g in Singular.independent_sets(Is)[1] for g in gens(Rs)]

    coeff_indets = typeof(R[1])[]
    remaining_indets = typeof(R[1])[]
    for i in 1:length(indep_indices)
        indep_indices[i] ? push!(coeff_indets, R[i]) : push!(remaining_indets, R[i])
    end
    gb = standard_basis(I, ordering=invlex([remaining_indets; coeff_indets]))

    S, vs = polynomial_ring(base_ring(R), Symbol.(coeff_indets))
    for (i,x) in enumerate(coeff_indets)
        :($(Symbol(x)) = $(vs[i])) |> eval
    end
    coeff_ring = fraction_field(S)

    Q, vq = polynomial_ring(coeff_ring, Symbol.(remaining_indets))
    for (i,x) in enumerate(remaining_indets)
        :($(Symbol(x)) = $(vq[i])) |> eval
    end
    I_in_Q = ideal([eval(Meta.parse("$(g)")) for g in gens(gb)])
    leading_coeffs = [leading_coefficient(g) for g in gens(I_in_Q)]
    [:($(Symbol(x)) = $(x)) |> eval for x in coeff_indets]
    h = Meta.parse("$(lcm(leading_coeffs))") |> eval
    Isat, sat_index = saturation_with_index(I, ideal(h))
    J = (I + ideal(h^sat_index)) : ideal(prod(elems))
    if is_one(J) || Isat ⊆ J
        return [I_in_Q]
    else
        return [I_in_Q; reduction_to_zero_dim(J, elems)]
    end
end

function exponent_lattice(I, elems)
    @assert I : ideal(elems) == I "`I` is not saturated wrt `elems`."
    @assert ~is_one(I) "`I` should be a proper ideal, but `I = (1)`."
    if dim(I) == 0
        return exponent_lattice_zero_dim(I, elems)
    else
        zero_dim_ideals = reduction_to_zero_dim(I, elems)
        lattice = integer_lattice(identity_matrix(ZZ, length(elems)))
        for J in zero_dim_ideals
            R = base_ring(J)
            K = coefficient_ring(J)
            for x in [gens(R); gens(K)]
               :($(Symbol(x)) = $(x)) |> eval
            end
            elems_in_J = R.([eval(Meta.parse("$(x)")) for x in elems])
            lattice = intersect(lattice, exponent_lattice_zero_dim(J, elems_in_J))
        end
        return lattice
    end
end


