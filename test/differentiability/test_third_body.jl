#! DIFFRACTOR FAILING AT TIME CONVERSION -- EXTRAPOLATION (SatelliteToolboxTransformations.jl)
#! ENZYME FAILING AT READ-ONLY ARG? (AstroForceModels.jl)
#! ZYGOTE FAILING A TUPLE SUBTRACTION (AstroForceModels.jl)

@testset "Third Body Differentiability State" begin

    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = ComponentVector(; JD=JD)
    t = 0.0

    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    for backend ∈ _BACKENDS
        if backend[1] == "Diffractor" || backend[1] == "Enzyme" || backend[1] == "Zygote"
            continue
        end
        testname = "Third Body Differentiability " * backend[1] 
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, p, t, sun_third_body),
                AutoFiniteDiff(),
                state
            )
        
            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, p, t, sun_third_body)),
                backend[2],
                state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad rtol=2e-1
        end
    end
end

@testset "Third Body Differentiability Time" begin

    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    eop_data = fetch_iers_eop()
    p = ComponentVector(; JD=JD)
    t = 0.0

    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    for backend ∈ _BACKENDS
        if backend[1] == "Diffractor" || backend[1] == "Enzyme" || backend[1] == "Zygote"
            continue
        end
        testname = "Third Body Differentiability " * backend[1] 
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(state, p, x, sun_third_body),
                AutoFiniteDiff(),
                t
            )
        
            f_ad, df_ad = value_and_derivative(
                (x) -> Array(acceleration(state, p, x, sun_third_body)),
                backend[2],
                t
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol=1e-10
        end
    end
end