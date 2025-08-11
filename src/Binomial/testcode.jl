using Oscar, Test
include("Binomial.jl")

@testset "BB Code" begin
    R, (x,y,s,t) = polynomial_ring(GF(2), [:x,:y,:s,:t])
    I = ideal([x^3+y^2+y, t^3+x^2+x, x*s-1, y*t-1])
end

@testset "X-Cube Code" begin
    R, (x,y,z,s,t,w) = polynomial_ring(GF(2), [:x,:y,:z,:s,:t,:w])
    I = ideal([(1+x)*(1+y),(1+y)*(1+z),(1+x)*(1+z),x*s-1,y*t-1,z*w-1])
end

@testset "Extension field" begin
    S, (x,) = polynomial_ring(GF(2), [:x])
    R, (y,z,s,t,w) = polynomial_ring(fraction_field(S), [:y,:z,:s,:t,:w])
    I = ideal([y*z + y + z + 1,
    y*s + y + s + 1,
    z*s + z + s + 1,
    y*t + 1,
    z*t + z + t + 1,
    s*t + s + t + 1,
    y*w + y + w + 1,
    z*w + 1,
    s*w + s + w + 1,
    t*w + t + w + 1,
    (x + 1)*y + x + 1,
    (x + 1)*z + x + 1,
    x*s + 1,
    (x + 1)*t + x + 1,
    (x + 1)*w + x + 1])
end