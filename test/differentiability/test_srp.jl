#! DIFFRACTOR FAILING AT TIME CONVERSION -- EXTRAPOLATION (SatelliteToolboxTransformations.jl)
#! ENZYME FAILING AT READ-ONLY ARG? (AstroForceModels.jl)
#! ZYGOTE FAILING WITH SOME KIND OF MISSHAPE? (AstroForceModels.jl)

@testset "SRP Differentiability State" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)
    t = 0.0

    eop_data = fetch_iers_eop()

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    satellite_srp_model = CannonballFixedSRP(0.2)

    sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    srp_model = SRPAstroModel(satellite_srp_model=satellite_srp_model, sun_data=sun_model, eop_data=eop_data, shadow_model=Conical())

    for backend ∈ _BACKENDS
        if backend[1] == "Diffractor" || backend[1] == "Enzyme" || backend[1] == "Zygote"
            continue
        end
        testname = "SRP Differentiability " * backend[1] 
        @testset "$testname" begin
    
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, p, t, srp_model),
                AutoFiniteDiff(),
                state
            )
        
            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, p, t, srp_model)),
                backend[2],
                state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol=1e-10
        end
    end
end


@testset "SRP Differentiability Time" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    eop_data = fetch_iers_eop()

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    satellite_srp_model = CannonballFixedSRP(0.2)

    sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    srp_model = SRPAstroModel(satellite_srp_model=satellite_srp_model, sun_data=sun_model, eop_data=eop_data, shadow_model=Conical())

    for backend ∈ _BACKENDS
        if backend[1] == "Diffractor" || backend[1] == "Enzyme" || backend[1] == "Zygote"
            continue
        end
        testname = "SRP Differentiability " * backend[1] 
        @testset "$testname" begin
    
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(state, p, x, srp_model),
                AutoFiniteDiff(),
                t
            )
        
            f_ad, df_ad = value_and_derivative(
                (x) -> Array(acceleration(state, p, x, srp_model)),
                backend[2],
                t
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol=1e-10
        end
    end
end


@testset "SRP Differentiability SRP Parameters" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    eop_data = fetch_iers_eop()

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    RC = .2

    sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

    for backend ∈ _BACKENDS
        if backend[1] == "Diffractor" || backend[1] == "Enzyme" || backend[1] == "Zygote"
            continue
        end
        testname = "SRP Differentiability " * backend[1] 
        @testset "$testname" begin
    
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(state, p, t, SRPAstroModel(satellite_srp_model=CannonballFixedSRP(x), sun_data=sun_model, eop_data=eop_data, shadow_model=Conical())),
                AutoFiniteDiff(),
                RC
            )
        
            f_ad, df_ad = value_and_derivative(
                (x) -> Array(acceleration(state, p, t, SRPAstroModel(satellite_srp_model=CannonballFixedSRP(x), sun_data=sun_model, eop_data=eop_data, shadow_model=Conical()))),
                backend[2],
                RC
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol=1e-10
        end
    end
end