using AllocCheck
using Aqua
using AstroForceModels
using ComponentArrays
using JET
using LinearAlgebra
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxCelestialBodies
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using SpaceIndices
using Test

using DifferentiationInterface
using FiniteDiff, ForwardDiff, Diffractor, Enzyme, Mooncake, PolyesterForwardDiff, Zygote

@testset "AstroForceModels.jl" begin
    # Drag Tests
    include("drag/test_satellite_shape_model.jl")
    include("drag/test_density_calculator.jl")
    include("drag/test_drag_accel.jl")

    # SRP Tests
    include("solar_radiation_pressure/test_satellite_shape_models.jl")
    include("solar_radiation_pressure/test_shadow_models.jl")
    include("solar_radiation_pressure/test_srp_accel.jl")

    # Third Body Tests
    include("third_body/test_celestial_body.jl")
    include("third_body/test_third_body_model.jl")
    include("third_body/test_third_body_accel.jl")

    # Zonal Harmonics Tests
    include("gravity/test_grav_accel.jl")

    # Relativity Tests
    include("relativity/test_relativity.jl")

    # Dynamics Builder
    include("test_dynamics_builder.jl")
end

using DifferentiationInterface, ForwardDiff, FiniteDiff, Mooncake

const _BACKENDS = (
    ("ForwardDiff", AutoForwardDiff()),
    ("Diffractor", AutoDiffractor()),
    ("Enzyme", AutoEnzyme()),
    ("Mooncake", AutoMooncake(; config=nothing)),
    ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
    ("Zygote", AutoZygote()),
)

@testset "Differentiability" begin
    include("differentiability/test_drag.jl")
    include("differentiability/test_srp.jl")
    include("differentiability/test_gravity.jl")
    include("differentiability/test_relativity.jl")
    include("differentiability/test_third_body.jl")
    include("differentiability/test_dynamics_builder.jl")
end

@testset "Aqua.jl" begin
    Aqua.test_all(AstroForceModels; ambiguities=(recursive = false))
end

@testset "JET Testing" begin
    rep = JET.test_package(
        AstroForceModels; toplevel_logger=nothing, target_modules=(@__MODULE__,)
    )
end

#TODO: GET THESE WORKING
#@testset "AllocCheck.jl" begin
#    # Force Model Allocation Check
#    include("test_allocations.jl")
#end
