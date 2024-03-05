# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Satellite Drag Models to easily compute the ballistic coefficient
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
Abstract satellite drag model to help compute drag forces.
"""
abstract type AbstractSatelliteDragModel end

"""
Cannonball Fixed Drag Model struct
Contains information to compute the ballistic coefficient of a cannonball drag model with a 
fixed ballistic coefficient.

# Fields
- `radius::Number`: The radius of the spacecraft.
- `mass::Number`: The mass of the spacecraft.
- `drag_coeff::Number`: The drag coefficient of the spacecraft.
- `ballistic_coeff::Number`: The fixed ballistic coefficient to use.
"""
struct CannonballFixedDrag <: AbstractSatelliteDragModel
    radius::Number
    mass::Number
    drag_coeff::Number
    ballistic_coeff::Number
end

"""
    CannonballFixedDrag(ballistic_coeff::Number)

Constructor for a fixed ballistic coefficient drag model.

# Arguments
- `ballistic_coeff::Number`: The fixed ballistic coefficient to use.

# Returns 
- `drag_model::CannonballFixedDrag`: A fixed ballistic coefficient drag model.

"""
function CannonballFixedDrag(ballistic_coeff::Number)
    return CannonballFixedDrag(1.0, 1.0, ballistic_coeff, ballistic_coeff)
end

"""
    CannonballFixedDrag(radius::Number, mass::Number, drag_coeff::Number)

Constructor for a fixed cannonball ballistic coefficient drag model.

The ballistic coefficient is computed with the following equation:

                BC = CD * area/mass

where area is the 2D projection of a sphere

                area = π * r^2

# Arguments
- `radius::Number`: The radius of the spacecraft.
- `mass::Number`: The mass of the spacecraft.
- `drag_coeff::Number`: The drag coefficient of the spacecraft.

# Returns 
- `drag_model::CannonballFixedDrag`: A fixed ballistic coefficient drag model.

"""
function CannonballFixedDrag(radius::Number, mass::Number, drag_coeff::Number)
    area = π * radius^2.0

    return CannonballFixedDrag(radius, mass, drag_coeff, drag_coeff * area / mass)
end

"""
    ballistic__coefficient(
        u::AbstractArray, 
        p::ComponentVector, 
        t::Number, 
        model::CannonballFixedDrag)

Returns the ballistic coeffient for a drag model given the model and current state of the simulation.

# Arguments
    - `u::AbstractArray`: The current state of the simulation.
    - `p::ComponentVector`: The parameters of the simulation.
    - `t::Number`: The current time of the simulation.
    - `model::CannonballFixedDrag`: Drag model for the spacecraft.

# Returns
    -`ballistic_coeff::Number`: The current ballistic coefficient of the spacecraft.

"""
@inline function ballistic_coefficient(
    u::AbstractArray, p::ComponentVector, t::Number, model::CannonballFixedDrag
)
    return model.ballistic_coeff
end

"""
State Drag Model struct
Empty struct used for when the simulation state includes the ballistic coefficient.
"""
struct StateDragModel <: AbstractSatelliteDragModel end

"""
    ballistic__coefficient(
        u::AbstractArray, 
        p::ComponentVector, 
        t::Number, 
        model::StateDragModel)

Returns the ballistic coeffient for a drag model given the model and current state 
of the simulation.

# Arguments
    - `u::AbstractArray`: The current state of the simulation.
    - `p::ComponentVector`: The parameters of the simulation.
    - `t::Number`: The current time of the simulation.
    - `model::StateDragModel`: The drag model for the spacecraft.

# Returns
    -`ballistic_coeff::Number`: The current ballistic coefficient of the spacecraft.
"""
@inline function ballistic_coefficient(
    u::AbstractArray, p::ComponentVector, t::Number, model::StateDragModel
)
    return u.BC
end
