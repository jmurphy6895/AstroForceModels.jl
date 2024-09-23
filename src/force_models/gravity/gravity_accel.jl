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
export GravityHarmonicsAstroModel, KeplerianGravityAstroModel

abstract type AbstractGravityAstroModel <: AbstractPotentialBasedForce end

"""
Gravitational Harmonics Astro Model struct
Contains information to compute the acceleration of a Gravitational Harmonics Model acting on a spacecraft.

# Fields
- `gravity_model::AbstractGravityModel`: The gravitational potential model and coefficient data.
- `eop_data::Union{EopIau1980,EopIau2000A}`: The data compute the Earth's orientation.
- `order::Int`: The maximum order to compute the graviational potential to, a value of -1 compute the maximum order of the supplied model. (Default=-1)
- `degree::Int`: The maximum degree to compute the graviational potential to, a value of -1 compute the maximum degree of the supplied model. (Default=-1)
"""
@with_kw struct GravityHarmonicsAstroModel{GT,EoT,V} <: AbstractGravityAstroModel where {
    GT<:AbstractGravityModel{<:Number},EoT<:Union{EopIau1980,EopIau2000A},V<:Int
}
    gravity_model::GT
    eop_data::EoT
    order::V = -1
    degree::V = -1
end

"""
    acceleration(u::AbstractArray, p::ComponentVector, t::Number, grav_model::GravityHarmonicsAstroModel)

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
    date = jd_to_date(p.JD + t / 86400.0)
    sec = floor(date[6])
    ms = floor((date[6] - floor(date[6])) * 1000.0)
    time = DateTime(date[1], date[2], date[3], date[4], date[5], sec, ms)
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

"""
Gravitational Keplerian Astro Model struct
Contains information to compute the acceleration of a Gravitational Harmonics Model acting on a spacecraft.

# Fields
- `μ::Number`: The gravitational potential constant of the central body.
"""
@with_kw struct KeplerianGravityAstroModel{MT} <:
                AbstractGravityAstroModel where {MT<:Number}
    μ::MT = μ_EARTH
end

"""
    acceleration(u::AbstractArray, p::ComponentVector, t::Number, grav_model::KeplerianGravityAstroModel)

Computes the gravitational acceleration acting on a spacecraft given a gravity model and current state and 
parameters of an object.

# Arguments
- `u::AbstractArray`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `gravity_model::KeplerianGravityAstroModel`: Gravity model struct containing the relevant information to compute the acceleration.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional gravity acceleration acting on the spacecraft.

"""
function acceleration(
    u::AbstractArray, p::ComponentVector, t::Number, grav_model::KeplerianGravityAstroModel
)
    r = SVector{3}(u[1:3])
    r_norm = norm(r)

    return SVector{3}((-grav_model.μ / (r_norm^3)) * r)
end
