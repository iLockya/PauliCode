
"""
    integral_closure(I)

Compute the integral closure of the ring `R/I`.
"""
function integral_closure(I)
    R = base_ring(I)
    A, _ = quo(R, I)
    (A̅, _, generators) = normalization(A, algorithm=:isPrime)[1]
    RA = base_ring(A̅)
    for k in 1:nvars(R)
        :($(Symbol(R[k])) = $(RA[k])) |> eval
    end
    module_gens = [generators[1], gens(generators[2])[1:end-1]...]

    return A̅.I, [RA(eval(Meta.parse("$(g)"))) for g in module_gens]
end

function compute_presentation(I, elems)
    R = base_ring(I) 
    fraction_field = R.base_ring # K(a1,...as)
    underlying_ring = fraction_field.base_ring # K[a1,...as]
    perfect_field = underlying_ring.base_ring # K

    var_name_y = ["y$k" for k in 1:length(elems)+1]
    var_name_underlying = [string(underlying_ring[k]) for k in 1:nvars(underlying_ring)]
    var_name_R = [string(R[k]) for k in 1:nvars(R)]
    var_names = cat(var_name_y, var_name_underlying, var_name_R; dims=1)
    R_elim, vs = polynomial_ring(perfect_field, var_names)
    for (i, x) in enumerate(var_names)
        :($(Symbol(x)) = $(vs[i])) |> eval
    end

    generators = gens(I)
    make_invertible = []
    for i in eachindex(generators)
        coeffs = coefficients(generators[i])
        denom_lcm = lcm([denominator(c) for c in coeffs])
        generators[i] = denom_lcm * generators[i]
        push!(make_invertible, eval(Meta.parse("$denom_lcm")))
    end
    for gen in generators
        push!(make_invertible, eval(Meta.parse("$gen")) |> leading_coefficient)
    end

    J = ideal([eval(Meta.parse(string(f))) for f in generators])
    k = length(elems)
    indeterminates = gens(R_elim)
    # add h_i*y_i-g_i*y_k+1 to J where f_j = g_i/h_i
    for i in eachindex(elems)
        coeffs = coefficients(elems[i])
        denom_lcm = lcm([denominator(c) for c in coeffs])
        cleared_denoms = denom_lcm * elems[i]
        denom_lcm = eval(Meta.parse(string(denom_lcm)))
        cleared_denoms = eval(Meta.parse(string(cleared_denoms)))
        J += ideal(denom_lcm * indeterminates[i] - cleared_denoms * indeterminates[k+1])
        push!(make_invertible, denom_lcm)
    end

    sq_free = one(R_elim)
    for (p, e) in factor(prod(make_invertible)).fac
        sq_free *= p
    end
    R_elim, J, _ = localize(J, [sq_free])
    indeterminates = gens(R_elim)
    Jelim = eliminate(J, indeterminates[k+2:end])
    Ry, vs = polynomial_ring(perfect_field, var_name_y)
    for (i, x) in enumerate(var_name_y)
        :($(Symbol(x)) = $(vs[i])) |> eval
    end
    return ideal(Ry, [eval(Meta.parse("$j")) for j in gens(Jelim)])
end

function exponent_lattice_extension_field(I, elems)
    @assert is_prime(I) "I is not prime."
    S_ideal = compute_presentation(I, elems)
    coeff_ring = coefficient_ring(S_ideal)
    S_closure_ideal, module_gens = integral_closure(S_ideal)
    S_closure_ring = base_ring(S_closure_ideal)
    t = gens(S_closure_ring)[length(elems)]
    a = gens(S_closure_ring)[1:length(elems)+1]
    primes = Set{MPolyIdeal}()
    for a_i in a
        union!(primes,Set(minimal_primes(S_closure_ideal + ideal(S_closure_ring, a_i))))
    end
    valuations = zeros(ZZRingElem, length(primes), length(elems))
    for (i, p) in enumerate(primes)
        vals = Int[]
        for f in a 
            r = 1
            while !((p^r+S_closure_ideal):(ideal(f)+S_closure_ideal) ⊆ p)
                r += 1
            end
            push!(vals, r-1)
        end
        valuations[i,:] = [vals[i] - vals[end] for i in 1:length(vals)-1]
    end
    ker = kernel(matrix(valuations), side=:right)
    if ker.c == 0
        return integer_lattice(ker)
    end
    fractions = []
    for k in 1:ker.c
        gen = ker[:,k]
        num = one(S_closure_ring)
        den = one(S_closure_ring)
        for i in eachindex(gen)
            if gen[i] < 0
                den *= S_closure_ring[i]^(-gen[i])
            else
                num *= S_closure_ring[i]^(gen[i])
            end
        end
        diff = total_degree(num) - total_degree(den)
        if diff > 0
            den *= t^diff
        elseif diff < 0
            num *= t^(-diff)
        end
        push!(fractions, [num, den])
    end
    field_elems = []
    for frac in fractions
        f = frac[1] * module_gens[1]
        repr_ideal = ideal([frac[1]*el for el in module_gens])+S_closure_ideal
        _, rep_coeffs = normal_form_with_transformation(f, repr_ideal)
        rep_coeffs = [c for c in rep_coeffs[1:length(module_gens)]]
        push!(field_elems, dot(rep_coeffs,[one(S_closure_ring),gens(S_closure_ring)[length(elems)+2:end]...]))
    end
    indeterminates = [string.(gens(S_closure_ring)); ["f$s" for s in 1:length(field_elems)]]
    elim_ring, vs = polynomial_ring(coeff_ring, indeterminates)
    for (i, x) in enumerate(vs)
        :($(Symbol(x)) = $(vs[i])) |> eval
    end
    field_ideal = ideal([eval(Meta.parse("$g")) for g in gens(S_closure_ideal)])
    field_ideal+= ideal([eval(Meta.parse("$(field_elems[end-i])"))-gens(elim_ring)[end-i] for i in 0:length(field_elems)-1])
    field_ideal = eliminate(field_ideal,vs[1:end-length(field_elems)])
    field_ring, vs = polynomial_ring(coeff_ring, ["f$s" for s in 1:length(field_elems)])
    for (i, x) in enumerate(vs)
        :($(Symbol(x)) = $(vs[i])) |> eval
    end
    field_ideal = ideal([eval(Meta.parse("$g")) for g in gens(field_ideal)])
    field_elems = [eval(Meta.parse("$g")) for g in last(gens(elim_ring), length(field_elems))]
    if characteristic(coeff_ring) == 0
        finite_extension_lattice = exponent_lattice_number_field_max(field_ideal, field_elems)
    else
        finite_extension_lattice = exponent_lattice_finite_field_max(field_ideal, field_elems)
    end
    lattice_gens = zeros(ZZRingElem, finite_extension_lattice.basis_matrix.r, length(elems))
    for k in 1:finite_extension_lattice.basis_matrix.r
        gen = finite_extension_lattice.basis_matrix[k,:]
        new_gen = sum([ZZ(gen[i]) * ker[:,i] for i in eachindex(gen)])   
        lattice_gens[k,:] = new_gen
    end
    return integer_lattice(matrix(lattice_gens))
end