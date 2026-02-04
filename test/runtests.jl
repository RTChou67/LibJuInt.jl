using LibJuInt
using Test

@testset "LibJuInt Integrals" begin
    # Define a simple Gaussian primitive (unnormalized for simplicity)
    # alpha = 0.5, coeff = 1.0, norm = 1.0
    pgtf = PGTF(0.5, 1.0, 1.0)
    
    # Define two s-type basis functions (L=0,0,0) at different positions
    # Basis 1 at Origin
    basis1 = Basis((0,0,0), [pgtf], (0.0, 0.0, 0.0))
    # Basis 2 at (1.0, 0.0, 0.0)
    basis2 = Basis((0,0,0), [pgtf], (1.0, 0.0, 0.0))

    # Define Atoms for Nuclear Attraction
    atom1 = Atom("H", 1, "Minimal", (0.0, 0.0, 0.0))
    atom2 = Atom("H", 1, "Minimal", (1.0, 0.0, 0.0))
    atoms = [atom1, atom2]

    @testset "Overlap Integral (S)" begin
        # Self overlap (unnormalized)
        # For alpha=0.5, p=1.0. Integral = (pi/p)^1.5 * exp(0) = pi^1.5 ≈ 5.568
        s11 = Sij(basis1, basis1)
        @test s11 isa Float64
        @test s11 > 0.0
        @test isapprox(s11, 5.56832799683170; atol=1e-5)

        # Overlap between centers
        s12 = Sij(basis1, basis2)
        @test s12 isa Float64
        @test s12 < s11 # Overlap should decrease with distance
        @test s12 > 0.0
    end

    @testset "Kinetic Energy Integral (T)" begin
        t11 = Tij(basis1, basis1)
        t12 = Tij(basis1, basis2)
        @test t11 isa Float64
        @test t12 isa Float64
        # Kinetic energy should be positive for self-interaction
        @test t11 > 0.0
    end

    @testset "Nuclear Attraction Integral (V)" begin
        v11 = Vij(basis1, basis1, atoms)
        v12 = Vij(basis1, basis2, atoms)
        @test v11 isa Float64
        @test v12 isa Float64
        # Attractive potential is negative
        @test v11 < 0.0 
    end

    @testset "Electron Repulsion Integral (G)" begin
        # (11|11)
        g1111 = Gijkl(basis1, basis1, basis1, basis1)
        @test g1111 isa Float64
        @test g1111 > 0.0

        # (11|22)
        g1122 = Gijkl(basis1, basis1, basis2, basis2)
        @test g1122 isa Float64
        @test g1122 > 0.0
        
        # Coulomb repulsion should decrease with distance (11|11) > (11|22)
        @test g1111 > g1122
    end
end
