
@testset "SRP Shadow Models" verbose = true begin
    eop_data = fetch_iers_eop()
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    R_MOD2J2000 = r_eci_to_eci(MOD(), J2000(), JD, eop_data)

    sun_pos = R_MOD2J2000 * sun_position_mod(JD)

    sat_pos = 6800 * normalize(sun_pos)

    shadow_model_set = [Cylindrical(), Conical(), Conical_Simplified(), INPE(), None()]

    # No Occulsion
    for model in shadow_model_set
        shadow_scaling = shadow_model(sat_pos, sun_pos, model; R_Sun=R_SUN, R_Earth=R_EARTH)
        @test shadow_scaling == 1.0
    end

    # Full Occulsion
    for model in shadow_model_set[1:(end - 1)]
        shadow_scaling = shadow_model(
            -sat_pos, sun_pos, model; R_Sun=R_SUN, R_Earth=R_EARTH
        )
        @test shadow_scaling == 0.0
    end
    @test shadow_model(-sat_pos, sun_pos, None(); R_Sun=R_SUN, R_Earth=R_EARTH) == 1.0

    # Partial Occulsion
    sun_pos_simple = [norm(sun_pos), 0.0, 0.0]
    sat_pos_simple = [-20000.0, 6378.1, 0.0]
    @test shadow_model(
        sat_pos_simple, sun_pos_simple, Cylindrical(); R_Sun=R_SUN, R_Earth=R_EARTH
    ) ≈ 0.0 atol=1E-9
    @test shadow_model(
        sat_pos_simple, sun_pos_simple, Conical(); R_Sun=R_SUN, R_Earth=R_EARTH
    ) ≈ 0.2573669402416321 atol=1E-9
    @test shadow_model(
        sat_pos_simple, sun_pos_simple, Conical_Simplified(); R_Sun=R_SUN, R_Earth=R_EARTH
    ) ≈ 0.303477350846976 atol=1E-9
    @test shadow_model(
        sat_pos_simple, sun_pos_simple, INPE(); R_Sun=R_SUN, R_Earth=R_EARTH
    ) ≈ 0.6874494027008462 atol=1E-9
    @test shadow_model(
        sat_pos_simple, sun_pos_simple, None(); R_Sun=R_SUN, R_Earth=R_EARTH
    ) ≈ 1.0 atol=1E-9
end
