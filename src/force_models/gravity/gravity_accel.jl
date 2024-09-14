# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Zonal Harmonics
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export GravityHarmonicsAstroModel

"""
SRP Astro Model struct
Contains information to compute the acceleration of a SRP a spacecraft.

# Fields
- `satellite_srp_model::AbstractSatelliteDragModel`: The satellite srp model for computing the ballistic coefficient.
- `sun_data::ThirdBodyModel`: The data to compute the Sun's position.
"""
@with_kw struct GravityHarmonicsAstroModel{GT,EoT,V} <: AbstractPotentialBasedForce where {
    GT<:AbstractGravityModel{<:Number},EoT<:Union{EopIau1980,EopIau2000A},V<:Int
}
    gravity_model::GT
    eop_data::EoT
    order::V = -1
    degree::V = -1
end

"""
    acceleration(u::AbstractArray, p::ComponentVector, t::Number, srp_model::GravityHarmonicsAstroModel)

Computes the gravitational acceleration acting on a spacecraft given a gravity model and current state and 
parameters of an object.

# Arguments
- `u::AbstractArray`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `gravity_model::GravityHarmonicsAstroModel`: Gravity model struct containing the relevant information to compute the acceleration.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional gravity acceleration acting on the spacecraft.

"""
function acceleration(
    u::AbstractArray, p::ComponentVector, t::Number, grav_model::GravityHarmonicsAstroModel
)
    # Compute the J2000 to ITRF rotation matrix
    R_J2002ITRF = r_eci_to_ecef(J2000(), ITRF(), p.JD + t / 86400.0, grav_model.eop_data)

    # Compute the ITRF position
    itrf_pos = R_J2002ITRF * @view(u[1:3]) .* 1E3

    #TODO: INTERNAL DATE REPRESENTATION TO JD
    # Compute the ITRF acceleration of the spacecraft
    time = DateTime(jd_to_date(p.JD + t / 86400.0)...)
    accel_itrf =
        GravityModels.gravitational_acceleration(
            grav_model.gravity_model,
            itrf_pos,
            time;
            max_degree=grav_model.degree,
            max_order=grav_model.order,
        ) ./ 1E3

    # Rotate into the J200 frame
    return R_J2002ITRF' * accel_itrf
end
