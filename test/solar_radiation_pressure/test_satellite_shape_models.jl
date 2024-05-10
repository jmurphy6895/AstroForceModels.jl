@testset "Cannonball Fixed SRP Satellite Shape" verbose = true begin
    radius = 10.0 * abs(randn())
    mass = 260.0 * abs(randn())
    reflectivity_coeff = 2.0 * abs(randn())

    area = π * (radius^2.0)
    RC = reflectivity_coeff * area / mass

    srp_model1 = CannonballFixedSRP(RC)
    srp_model2 = CannonballFixedSRP(radius, mass, reflectivity_coeff)

    u = randn(6)
    p = [0.0]
    t = 1000.0 * abs(randn())

    @test reflectivity_ballistic_coefficient(u, p, t, srp_model1) == RC
    @test reflectivity_ballistic_coefficient(u, p, t, srp_model2) == RC

    @test srp_model1.radius == 1.0
    @test srp_model1.mass == 1.0
    @test srp_model1.reflectivity_coeff == RC

    @test srp_model2.radius == radius
    @test srp_model2.mass == mass
    @test srp_model2.reflectivity_coeff == reflectivity_coeff
end

@testset "Cannonball Fixed SRP Satellite Shape Errors" verbose = true begin
    @test_throws ArgumentError CannonballFixedSRP(-1.0)
    @test_throws ArgumentError CannonballFixedSRP(-10.0, 260.0, 1.0)
    @test_throws ArgumentError CannonballFixedSRP(10.0, -260.0, 1.0)
    @test_throws ArgumentError CannonballFixedSRP(10.0, 260.0, -1.0)
end
