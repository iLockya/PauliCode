
"""
    factor_multiplicity(n, p)

Find the largest integer `k` such that `p^k` divides `n`.
"""
function factor_multiplicity(n,p)
    for (f, e) in factor(n).fac
        f == p && return e
    end
    return 0
end

function separable_part(f, levels)
    h = gcd(f, derivative(f))
    g1 = f / h
    newH = gcd(h, derivative(h))
    while newH != h
        h = newH
        newH = gcd(h, derivative(h))
    end
    if h == 1
        return g1, levels
    else
        p = characteristic(f.parent)
        k = 0
        for (k, ci) in h.coeffs
            ci == 0 && continue
            ei = k - 1 # h = ci * a^ei
            if ei != 0
                factor_mul = factor_multiplicity(ei, p)
                if factor_mul == 0
                    k = 0
                    break
                elseif factor_mul != 0 && k == 0
                    k = factor_mul
                else
                    k = min(k, factor_mul)
                end
            end
        end
        g2 = one(f.parent)
        if k ≥ 1
            g2, levels = pkth_root(h, k, levels)
        end
        return separable_part(g1*g2, levels)
    end 
end

"""
    scale(f, scale_levels, levels)

# Arguments
- `f`: A polynomial in `GF(p)(t_1,...t_s)[x_1,...x_n]`.
"""
function scale(f, scale_levels, levels)
    R = f.parent
    scaled_f = zero(R)
    if R isa PolyRing
        for (i, coeff) in enumerate(f.coeffs)
            mon = gen(R)^(i-1)
            scaled_numerator = scale_polynomial(coeff.num, scale_levels, levels)
            scaled_denominator = scale_polynomial(coeff.den, scale_levels, levels)
            scaled_f += (scaled_numerator / scaled_denominator) * mon
        end
    elseif R isa MPolyRing
        for (mon, coeff) in zip(monomials(f), coefficients(f))
            scaled_numerator = scale_polynomial(coeff.num, scale_levels, levels)
            scaled_denominator = scale_polynomial(coeff.den, scale_levels, levels)
            scaled_f += (scaled_numerator / scaled_denominator) * mon
        end
    end
    return scaled_f
end

function scale_polynomial(f, scale_levels, levels)
    R = f.parent
    p = characteristic(R)
    scaled_f = zero(R)
    for (mon, coeff) in zip(monomials(f), coefficients(f))
        exponents_list = collect(exponents(mon))[1]
        for i in eachindex(exponents_list)
            if exponents_list[i] != 0 && levels[i] < scale_levels[i]
                exponents_list[i] = exponents_list[i]*p^(scale_levels[i]-levels[i])
            end
        end
        scaled_f += coeff * prod(gens(R).^exponents_list)
    end
    return scaled_f
end


function compute_decomposition(I, elems, sep_parts)
    R = base_ring(I)
    A, p = quo(R, I)
    decompositions = []

    for i in 1:length(elems)
        f = sep_parts[i]
        g = derivative(f)
        sep = elems[i] |> p
        nil = 0
        while f(sep) != 0
            nil += f(sep) * g(sep)^(-1)
            sep -= f(sep) * g(sep)^(-1)
        end
        push!(decompositions,lift.([A(sep), A(nil)]))
    end
    return decompositions
end

"""
    separable_decomposition_perfect(I, elems)

# Arguments
- `I`: A zero-dimensional ideal over a perfect ring.
- `elems`: A vector of elements in the ring.
"""
function separable_decomposition_perfect(I, elems)
    sep_parts = []
    for f in elems
        min_poly = minimal_polynomial(I, f)
        sq_free_part = factor(min_poly).fac |> keys |> prod
        push!(sep_parts, sq_free_part)
    end
    return compute_decomposition(I, elems, sep_parts)
end

function separable_decomposition_non_perfect(I, elems)
    coeff_base = base_ring(base_ring(I))
    unscaled_parts = []
    levels_list = []
    n = ngens(coeff_base)
    for f in elems
        μ = minimal_polynomial(I, f)
        levels = repeat([0], n)
        sep_part, levels = separable_part(μ, levels)
        push!(unscaled_parts, sep_part)
        push!(levels_list, levels)
    end
    max_levels = repeat([0], n)
    for i in 1:n
        max_levels[i] = maximum([levels[i] for levels in levels_list])
    end
    sep_parts = [scale(unscaled_parts[i], max_levels, levels_list[i]) for i in 1:length(elems)]
    scaled_gens = [scale(f, max_levels, repeat([0], n)) for f in gens(I)]
    elems = [scale(f, max_levels, repeat([0], n)) for f in elems]
    scaledI = ideal(scaled_gens)
    return compute_decomposition(scaledI, elems, sep_parts), scaledI
end
"""
    nil_index(I, f)

Return the minimul integer `k>0`, such that `f^k ∈ I`. 
# Examples
```
```
"""
function nil_index(I, f)
    k = 1
    fk = f
    while fk ∉ I # TODO: check f is nilponent.
        fk *= f
        k  += 1
    end
    k
end

function extended_discrete_log(f, order_f, gens, orders, I)
    for i in 1:order_f-1
        combination = zeros(Int, length(gens))
        while has_next_combination(combination, orders)
            combination = next_combination(combination, orders)
            group_elem = one(base_ring(I))
            for j in 1:length(gens)
                group_elem *= normal_form(gens[j]^combination[j], I)
            end
            f_i = normal_form(f^i, I)
            if f_i*group_elem-1 in I
                return [i;combination]
            end
        end
    end
    return [order_f; zeros(Int, length(gens))]
end

"""
    unipotent_lattice_p(I, elems)

Return the exponent lattice of `elems` in the unipotent subgroup `1+Rad(0)` of `(R/I)ˣ`.
"""
function unipotent_lattice_p(I, elems)
    p = characteristic(base_ring(I))
    orders = ZZRingElem[]
    for f in elems
        index = nil_index(I, f-1)
        logf = Int(ceil(log(index)/log(p)))
        push!(orders, p^logf) # ord(1+f) = p^k, where k = min{k: p^k ≥ ord(f)}
    end
    last_row = zeros(ZZRingElem, length(elems))
    last_row[end] = orders[end]
    lattice_rows = [last_row]
    for i in length(elems)-1:-1:1
        ext_log = extended_discrete_log(elems[i], orders[i], elems[i+1:end], orders[i+1:end], I)
        orders[i] = ext_log[1]
        push!(lattice_rows, [zeros(ZZRingElem, i-1); ext_log])
    end
    return integer_lattice(matrix(lattice_rows))
end


function exponent_lattice_zero_dim(I, elems)
    @assert dim(I) == 0
    R = base_ring(I)
    K = coefficient_ring(I)
    _, proj = quo(R, I)

    if characteristic(K) == 0 || K isa FqField
        decompositions = separable_decomposition_perfect(I, elems)
    else
        decompositions, I = separable_decomposition_non_perfect(I, elems)
    end

    nil_image = MPolyRingElem[]
    sep_parts = MPolyRingElem[]
    for i in 1:length(elems)
        sep = decompositions[i][1]
        push!(sep_parts, sep)
        push!(nil_image, lift(elems[i] * inv(proj(sep))))
    end
    if characteristic(K) == 0
        # TODO
        lattice = unipotent_lattice_0(I, nil_image)
    else
        lattice = unipotent_lattice_p(I, nil_image)
    end
    #is_empty(lattice) && return lattice

    for max_ideal in minimal_primes(I)
        if K == QQ
            # TODO
            exp_lattice = exponent_lattice_number_field_max(max_ideal, sep_parts)
        elseif K isa FqField
            exp_lattice = exponent_lattice_finite_field_max(max_ideal, sep_parts)
        else
            exp_lattice = exponent_lattice_extension_field(max_ideal, sep_parts)
        end
        #is_empty(exp_lattice) && return exp_lattice
        lattice = intersect(lattice, exp_lattice) # TODO
    end

    return lattice
end