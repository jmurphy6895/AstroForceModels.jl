@testset "Test Celestial Body Construction" begin
    
    sun = Sun()
    @test sun.μ == 1.32712440018E11
    @test sun.name == :Sun
    @test sun.Req == 6.955E5

    moon = Moon()
    @test moon.μ == 4.90486959E3
    @test moon.name == :Moon
    @test moon.Req == 1738.1

    sun_2 = CelestialBody(:Sun)
    moon_2 = CelestialBody(:Moon)

    @test sun.name == sun_2.name
    @test moon.name == moon_2.name

end

@testset "Test Low Precision Celestial Body Construction"  begin
    
    sun = Sun(; T=Float32)
    @test sun.μ == Float32(1.32712440018E11)
    @test sun.name == :Sun
    @test sun.Req == Float32(6.955E5)

    moon = Moon(; T=Float32)
    @test moon.μ == Float32(4.90486959E3)
    @test moon.name == :Moon
    @test moon.Req == Float32(1738.1)

end

@testset "Test Custom Celestial Body Construction" begin
    
    tatooine = CelestialBody(
        :Tatooine,
        :Dual_Sun,
        jpl_code=99999,
        μ=1.2E8,
        Req=3456,
    )

    @test tatooine.name == :Tatooine
    @test tatooine.central_body ==:Dual_Sun
    @test tatooine.jpl_code == 99999
    @test μ == 1.2E8
    @test Req == 3456

end

@testset "Test Celestial Body Errors" begin
    
    @test_throws ArgumentError CelestialBody(:Kerbal)
    
end