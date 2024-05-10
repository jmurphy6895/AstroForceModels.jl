using AstroForceModels
using ComponentArrays
using LinearAlgebra
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxCelestialBodies
using SatelliteToolboxTransformations
using SpaceIndices
using Test

@testset "AstroForceModels.jl" begin
    # Drag Tests
    include("drag/test_satellite_shape_model.jl")
    include("drag/test_density_calculator.jl")
    include("drag/test_drag_accel.jl")

    # SRP Tests
    include("solar_radiation_pressure/test_satellite_shape_models.jl")
    include("solar_radiation_pressure/test_shadow_models.jl")
    include("solar_radiation_pressure/test_srp_accel.jl")
end
