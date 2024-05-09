using AstroForceModels
using ComponentArrays
using LinearAlgebra
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxTransformations
using SpaceIndices
using Test

@testset "AstroForceModels.jl" begin
    # Drag Tests
    include("drag/test_satellite_shape_model.jl")
    include("drag/test_density_calculator.jl")
    include("drag/test_drag_accel.jl")

    # SRP Tests
end
