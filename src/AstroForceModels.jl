module AstroForceModels

using AngleBetweenVectors
using ComponentArrays, StaticArraysCore
using LinearAlgebra
using Parameters
using SatelliteToolboxBase
using SatelliteToolboxCelestialBodies
using SatelliteToolboxGravityModels
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxTransformations
using SpaceIndices
using ValSplit

abstract type AbstractAstroForceModel end

abstract type AbstractNonPotentialBasedForce <: AbstractAstroForceModel end
abstract type AbstractPotentialBasedForce <: AbstractAstroForceModel end

include("./constants.jl")

include("./force_models/drag/satellite_shape_model.jl")
include("./force_models/drag/density_calculator.jl")
include("./force_models/drag/drag_accel.jl")

include("./force_models/third_body/celestial_body.jl")
include("./force_models/third_body/third_body_model.jl")
include("./force_models/third_body/third_body_accel.jl")

include("./force_models/relativity/relativity_accel.jl")

include("./force_models/solar_radiation_pressure/satellite_shape_model.jl")
include("./force_models/solar_radiation_pressure/shadow_models.jl")
include("./force_models/solar_radiation_pressure/srp_accel.jl")

include("./force_models/gravity/gravity_accel.jl")

export acceleration

end
