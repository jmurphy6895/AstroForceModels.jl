#! DIFFRACTOR FAILING AT TIME CONVERSION -- EXTRAPOLATION (SatelliteToolboxTransformations.jl)
#! ENZYME FAILING AT READ-ONLY ARG? (AstroForceModels.jl)
#! ZYGOTE FAILING AT DCM CONSTRUCTOR (ReferenceFrameRotations.jl)
const _ATMOSPHERE_MODELS = (
    ("JB2008", JB2008()),
    ("JR1971", JR1971()),
    ("MSIS2000", MSIS2000()),
    ("ExpAtmo", ExpAtmo()),
    ("None", None()),
)

# See SatelliteToolboxAtmosphericModels.jl for details
const _SKIP_TESTS_DIFFRACTOR = ["JR1971", "MSIS2000", "JB2008"]
const _SKIP_TESTS_ENZYME = ["JR1971", "MSIS2000"]
const _SKIP_TESTS_ZYGOTE = ["JR1971", "MSIS2000"]

@testset "Drag Differentiability State" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)
    t = 0.0

    SpaceIndices.init()
    eop_data = fetch_iers_eop()

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    satellite_drag_model = CannonballFixedDrag(0.2)

    for backend in _BACKENDS
        if backend[1] == "Diffractor" || backend[1] == "Enzyme" || backend[1] == "Zygote"
            continue
        end
        testname = "Drag Differentiability " * backend[1]
        @testset "$testname" begin
            for atmo in _ATMOSPHERE_MODELS
                if (backend[1] == "Diffractor" && atmo[1] ∈ _SKIP_TESTS_DIFFRACTOR) ||
                    (backend[1] == "Enzyme" && atmo[1] ∈ _SKIP_TESTS_ENZYME) ||
                    (backend[1] == "Zygote" && atmo[1] ∈ _SKIP_TESTS_ZYGOTE)
                    continue
                end

                drag_model = DragAstroModel(satellite_drag_model, atmo[2], eop_data)

                f_fd, df_fd = value_and_jacobian(
                    (x) -> acceleration(x, p, t, drag_model), AutoFiniteDiff(), state
                )

                f_ad, df_ad = value_and_jacobian(
                    (x) -> Array(acceleration(x, p, t, drag_model)), backend[2], state
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad rtol = 2e-1
            end
        end
    end
    SpaceIndices.destroy()
end

@testset "Drag Differentiability Time" begin
    JD = date_to_jd(2024, 1, 5, 12, 32, 0.0)
    p = ComponentVector(; JD=JD)
    t = 0.0

    SpaceIndices.init()
    eop_data = fetch_iers_eop()

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    satellite_drag_model = CannonballFixedDrag(0.2)

    for backend in _BACKENDS
        if backend[1] == "Diffractor" || backend[1] == "Enzyme" || backend[1] == "Zygote"
            continue
        end
        testname = "Drag Differentiability " * backend[1]
        @testset "$testname" begin
            for atmo in _ATMOSPHERE_MODELS
                if (backend[1] == "Diffractor" && atmo[1] ∈ _SKIP_TESTS_DIFFRACTOR) ||
                    (backend[1] == "Enzyme" && atmo[1] ∈ _SKIP_TESTS_ENZYME) ||
                    (backend[1] == "Zygote" && atmo[1] ∈ _SKIP_TESTS_ZYGOTE)
                    continue
                end

                drag_model = DragAstroModel(satellite_drag_model, atmo[2], eop_data)

                f_fd, df_fd = value_and_derivative(
                    (x) -> acceleration(state, p, x, drag_model), AutoFiniteDiff(), t
                )

                f_ad, df_ad = value_and_derivative(
                    (x) -> Array(acceleration(state, p, x, drag_model)), backend[2], t
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad atol = 2e-1
            end
        end
    end
    SpaceIndices.destroy()
end

@testset "Drag Differentiability Model Parameters" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)
    t = 0.0

    SpaceIndices.init()
    eop_data = fetch_iers_eop()

    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    BC = 0.2

    for backend in _BACKENDS
        if backend[1] == "Diffractor" || backend[1] == "Enzyme" || backend[1] == "Zygote"
            continue
        end
        testname = "Drag Differentiability " * backend[1]
        @testset "$testname" begin
            for atmo in _ATMOSPHERE_MODELS
                if (backend[1] == "Diffractor" && atmo[1] ∈ _SKIP_TESTS_DIFFRACTOR) ||
                    (backend[1] == "Enzyme" && atmo[1] ∈ _SKIP_TESTS_ENZYME) ||
                    (backend[1] == "Zygote" && atmo[1] ∈ _SKIP_TESTS_ZYGOTE)
                    continue
                end

                f_fd, df_fd = value_and_derivative(
                    (x) -> acceleration(
                        state,
                        p,
                        t,
                        DragAstroModel(CannonballFixedDrag(x), atmo[2], eop_data),
                    ),
                    AutoFiniteDiff(),
                    BC,
                )

                f_ad, df_ad = value_and_derivative(
                    (x) -> Array(
                        acceleration(
                            state,
                            p,
                            t,
                            DragAstroModel(CannonballFixedDrag(x), atmo[2], eop_data),
                        ),
                    ),
                    backend[2],
                    BC,
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad rtol = 1e-4
            end
        end
    end
    SpaceIndices.destroy()
end
