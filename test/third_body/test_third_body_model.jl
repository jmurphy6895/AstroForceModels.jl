@testset "Third Body Model Ephemeris" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()

    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

    sun_pos = sun_third_body(JD, Position())

    expected_sun_pos = [36327254721.23421, -130787280761.2391, -56694897390.09085]

    @test expected_sun_pos ≈ sun_pos rtol = 1E-4

    sun_vel = sun_third_body(JD, Velocity())

    expected_sun_vel = [29356.8910247322, 6844.1181149359, 2965.5158876916]

    @test expected_sun_vel ≈ sun_vel rtol = 1E-2

    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    moon_pos = moon_third_body(JD, Position())

    expected_moon_pos = [-344989050.2810446, -175734543.4089319, -81943013.49241804]

    @test expected_moon_pos ≈ moon_pos rtol = 1E-3
end
