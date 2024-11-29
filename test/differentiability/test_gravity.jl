#! DIFFRACTOR FAILING AT TIME CONVERSION -- EXTRAPOLATION (SatelliteToolboxTransformations.jl)
#! ENZYME FAILING AT READ-ONLY ARG? (AstroForceModels.jl)
#! ZYGOTE FAILING A RESHAPE (AstroForceModels.jl)
@testset "Harmonics Differentiability State" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)
    t = 0.0

    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    for backend in _BACKENDS
        if backend[1] == "Diffractor" || backend[1] == "Enzyme" || backend[1] == "Zygote"
            continue
        end
        testname = "Gravity Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> acceleration(x, p, t, grav_model), AutoFiniteDiff(), state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(acceleration(x, p, t, grav_model)), backend[2], state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad rtol = 2e-1
        end
    end
end

@testset "Harmonics Differentiability Time" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)
    t = 0.0

    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    for backend in _BACKENDS
        if backend[1] == "Diffractor" || backend[1] == "Enzyme" || backend[1] == "Zygote"
            continue
        end
        testname = "Gravity Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> acceleration(state, p, x, grav_model), AutoFiniteDiff(), t
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(acceleration(state, p, x, grav_model)), backend[2], t
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-10
        end
    end
end
