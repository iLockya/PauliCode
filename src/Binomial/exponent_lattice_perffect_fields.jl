
"""
    primitive_element(I)
    primitive_element(I, B, group_order)

Return a primitive element (generator of multiplicative group) of the finite field `GF(p)[x_1,...x_n]/I`.
The primitive element is not unique.

# Arguments
- `I`: A maximum ideal in `GF(p)[x_1,...x_n]`;
- `B`: A basis of polynomials representing a basis of the vector space `GF(p)[x_1,...x_n]/I` over `GF(p)`;
- `group_order`: The order of the multiplicative group of `GF(p)[x_1,...x_n]/I`.

# Examples
```jldoctest
julia> S,(x,) = polynomial_ring(GF(3), [:x]);
julia> I = ideal(S, x^3+2x+1); # x^3+2x+1 is irreducible.
julia> primitive_element_finite(I)
2*x^2 + x
julia> primitive_element_finite(I)
x^2 + 2*x
julia> primitive_element_finite(I)
x^2 + 2
```
"""
function primitive_element(I)
    A, _ = quo(base_ring(I), I)
    p = characteristic(A)
    B = monomial_basis(A)
    group_order = p^length(B)-1
    return primitive_element(I, B, group_order)
end

function primitive_element(I, B, group_order)
    R = base_ring(I)
    p = characteristic(R)
    primitive_elem = one(R)
    for (fac, pow) in factor(group_order).fac
        found = false
        while !found
            random_elem = zero(R)
            for b in B
                random_elem += rand(0:p-1) * b
            end
            if random_elem ∉ I
                power = power_reduce(random_elem, group_order÷fac, I)
                found = (power-1) ∉ I
                found && (primitive_elem *= power_reduce(random_elem, group_order÷(fac^pow), I))
            end
        end
    end
    primitive_elem
end

"""
    repr_primitive_element(I, primitive_elem, elements, B)

Represent each element in `elements` with respect to the `primitive_elem` in `GF(p)[x_1,...x_n]/I`.

# Arguments
- `I`: A maximum ideal in `GF(p)[x_1,...x_n]`;
- `primitive_elem`: A primitive element in `GF(p)[x_1,...x_n]/I`;
- `elements`: A vector of elements in `I`;
- `B`: A basis of polynomials representing a basis of the vector space `GF(p)[x_1,...x_n]/I` over `GF(p)`.
"""
function repr_primitive_element(I, primitive_elem, elements, B)
    d = minimal_polynomial(I, primitive_elem) |> degree
    K = coefficient_ring(I)
    primitive_powers = zeros(K, d, length(B))
    for i in 1:d
        primitive_i = power_reduce(primitive_elem, i-1, I)
        D = Dict(monomials(primitive_i) .=> coefficients(primitive_i))
        for (j, b) in enumerate(B)
            primitive_powers[i,j] = haskey(D, b) ? K(D[b]) : zero(K)
        end
    end
    primitive_powers = matrix(primitive_powers)
    solutions = []
    for elem in elements
        el_repr = zeros(K, length(B))
        elem = normal_form(elem, I)
        D = Dict(monomials(elem) .=> coefficients(elem))
        for (j, b) in enumerate(B)
            el_repr[j] = haskey(D, b) ? K(D[b]) : zero(K)
        end
        push!(solutions, solve(primitive_powers, el_repr))
    end
    solutions 
end

"""
    exponent_lattice_finite_field(K, elements)

Return the exponent lattice `[k1,...km]` of `elements: [e1,...em]`. `e1^k1···en^kn=1`.

# Arguments
- `K`: A finite field;
- `elements`: `[e_1,...e_m]` in `K`.
"""
function exponent_lattice_finite_field(K, elements)
    group_order = order(K) - 1
    primitive_elem = primitive_element(K)
    logs = [disc_log(primitive_elem, elems) for elems in elements]
    push!(logs, -group_order)
    solution = kernel(matrix(logs))
    return integer_lattice(solution[:,1:end-1])
end

"""
    exponent_lattice_finite_field_max(I, elements)

Return the exponent lattice `[k1,...km]` of `elements: [e1,...em]`. `e1^k1···en^kn=1 mod I`.

# Arguments
- `I`: A maximum ideal in `GF(p)[x_1,...x_n]`;
- `elements`: `[e_1,...e_m]` in `I`.
"""
function exponent_lattice_finite_field_max(I, elements)
    R = base_ring(I)
    @assert is_finite(coefficient_ring(I)) "Coefficient field must be finite."
    p = characteristic(R)
    A, pmap = quo(R, I)
    B = monomial_basis(A)
    group_order = p^length(B)-1
    primitive_elem = primitive_element(I, B, group_order)
    μ = minimal_polynomial(I, primitive_elem)
    K = GF(μ) # field extension
    representations = repr_primitive_element(I, primitive_elem, elements, B)
    elems_wrt_gen = FqFieldElem[]
    for repr in representations
        repr_elem = [gen(K)^(i-1)*repr[i] for i in 1:length(repr)] |> sum
        push!(elems_wrt_gen, repr_elem)
    end
    return exponent_lattice_finite_field(K, elems_wrt_gen)
end

function exponent_lattice_number_field_max(I, elems)
end