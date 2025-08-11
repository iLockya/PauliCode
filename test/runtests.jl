using PauliCode
using Test

@testset "[144,16,16] BB code" begin
    C, (x,y) = CSSCode(["x^3+y^2+y","y^3+x^2+x"],[12,12])
    @test physical_qubits(C) == 144
    @test logical_qubits(C) == 16
    @test topological_index(C) == 8
end

@testset "BB" begin
    C, (x,y) = CSSCode(["x^-3+y^2+y","y^3+x^2+x"],[12,12])
end

@testset "Fish TC" begin
    C, (x,y) = CSSCode(["1+x^-1+x^-2*y+x^-2*y^-1","1+y^-1"])
end


@testset "X_Cube model" begin
    X_Cube, (x,y,z) = CSSCode(["1+x^-1+y^-1+(x*y)^-1" "1+y^-1+z^-1+(y*z)^-1" "1+x^-1+z^-1+(x*z)^-1"], 
                              ["1+z" "1+x" "0"; "0" "1+x" "1+y"], 
                              ["x^9-1","y^9-1","z^9-1"])
    @test GSD(X_Cube) == 51
end

@testset "polynomial_period" begin
    R, (x,y) = polynomial_ring(GF(2), [:x,:y])
    S, a = polynomial_ring(GF(3), :a)
    f = x^11+x^10+x^9+x^7+x^6+x^4+1
    g = a^2+a+2
    @test polynomial_period(f) == 762
    @test polynomial_period(g) == 8
end

@testset "Mobility sub lattice" begin
    C, (x,y) = CSSCode(["1+x+x^(-1)*y^(-3)","1+y+y^(-1)*x^3"],["x^12-1","y^12-1"]);
    @test MobilitySubLattice(C) == MobilityLattice([762, 762])
end
