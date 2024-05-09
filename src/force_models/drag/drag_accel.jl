# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Drag
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#TODO: REFERENCE
#   [1] 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export DragAstroModel, drag_accel

"""
Drag Astro Model struct
Contains information to compute the acceleration of a drag force on a spacecraft.

# Fields
- `satellite_drag_model::AbstractSatelliteDragModel`: The satellite drag model for computing the ballistic coefficient.
- `atmosphere_model::Symbol`: The atmospheric model for computing the density.
- `eop_data::EopIau1980`: Earth orientation parameters to help compute the density with the atmospheric model.
"""
struct DragAstroModel{T,V} <: AbstractNonPotentialBasedForce where {
    T<:AbstractSatelliteDragModel,V<:Union{EopIau1980,EopIau2000A}
}
    satellite_drag_model::T
    atmosphere_model::Symbol
    eop_data::V
end

"""'
    acceleration(drag_model::DragAstroModel, u::AbstractArray, p::ComponentVector, t::Number)

Computes the drag acceleration acting on a spacecraft given a drag model and current state and 
parameters of an object.

# Arguments
- `drag_model::DragAstroModel`: Drag model struct containing the relevant information to compute the acceleration.
- `u::AbstractArray`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional drag acceleration acting on the spacecraft.

"""
function acceleration(
    u::AbstractArray, p::ComponentVector, t::Number, drag_model::DragAstroModel
)
    # Compute density at the satellite's current position
    rho = compute_density(
        p.JD + t / 86400.0, u, drag_model.eop_data, drag_model.atmosphere_model
    )

    #TODO: OFFER OPTION TO COMPUTE FROM EOP or SPICE EPHEMERIS 
    ω_vec = SVector{3}([0.0; 0.0; EARTH_ANGULAR_SPEED])

    # Compute the ballistic coefficient
    BC = ballistic_coefficient(u, p, t, drag_model.satellite_drag_model)

    # Return the 3-Dimensional Drag Force
    return drag_accel(u, rho, BC, ω_vec, t)
end

"""
    drag_accel(u::AbstractArray, rho::Number, BC::Number, ω_vec::AbstractArray, t::Number, [DragModel]) -> SVector{3}{Number}

Compute the Acceleration Atmospheric Drag

The atmosphere is treated as a solid revolving with the Earth and the apparent velocity of the satellite is computed
using the transport theorem

                𝐯_app = 𝐯 - 𝛚 x 𝐫

The acceleration from drag is then computed with a cannonball model as
                
                𝐚 = 1/2 * ρ * BC * |𝐯_app|₂^2 * v̂


!!! note
    Currently only fixed cannonball state based ballistic coefficients are supported, custom models can be created for
    higher fidelity.

# Arguments

- `u::AbstractArray`: The current state of the spacecraft in the central body's inertial frame.
- `rho::Number`: Atmospheric density at (t, u) [kg/m^3].
- `BC::Number`: The ballistic coefficient of the satellite -- (area/mass) * drag coefficient [kg/m^2].
- `ω_vec::AbstractArray`: The angular velocity vector of Earth. Typically approximated as [0.0; 0.0; ω_Earth]
- `t::Number`: Current time of the simulation.

# Returns

- `SVector{3}{Number}`: Inertial acceleration from drag
"""
@inline function drag_accel(
    u::AbstractArray, rho::Number, BC::Number, ω_vec::AbstractArray, t::Number
)

    # Compute Apparent Velocity w.r.t the Atmosphere using the Transport Theorem
    apparent_vel = @view(u[4:6]) - cross(ω_vec, @view(u[1:3]))

    # Scaled by 1E3 to convert to km/s
    # TODO: HANDLE UNITS BETTER
    return SVector{3}((-0.5 * BC * rho * norm(apparent_vel) * apparent_vel) .* 1E3)
end
