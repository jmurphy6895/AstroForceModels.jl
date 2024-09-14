# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Third Body Model and Ephemeris Functions
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export Vallado
abstract type AbstractEphemerisType end
struct Vallado <: AbstractEphemerisType end

export Position, Velocity
abstract type EphemerisReturn end
struct Position <: EphemerisReturn end
struct Velocity <: EphemerisReturn end

export ThirdBodyModel
"""
Third Body Model Astro Model struct
Contains information to compute the acceleration of a third body force acting on a spacecraft.

# Fields
- `body::CelestialBody`: Celestial body acting on the craft.
- `ephem_type::AbstractEphemerisType`: Ephemeris type used to compute body's position. Options are currently Vallado().
"""
@with_kw struct ThirdBodyModel{BT,EoT,EpT} <: AbstractNonPotentialBasedForce where {
    BT<:CelestialBody,EoT<:Union{EopIau1980,EopIau2000A,Nothing},EpT<:AbstractEphemerisType
}
    body::BT = SunBody
    eop_data::EoT = nothing
    ephem_type::EpT = Vallado()
end

#TODO: EXPAND TO SPICE WITH EXTENSIONS
"""
Computes the position of the celestial body using Vallado's ephemeris

# Arguments
- `ephem_type::Symbol`: Ephemeris type used to compute body's position.
- `body::CelestialBody`: Celestial body acting on the craft.
- `time::Number`: Current time of the simulation in seconds.

# Returns
- `body_position: SVector{3}`: The 3-dimensional third body position in the J2000 frame.
"""
function get_position(
    ephem_type::Vallado, body::CelestialBody, eop_data::T, time::TT
) where {T<:Union{EopIau1980,EopIau2000A,Nothing},TT}

    # Compute the MOD frame in the J2000 frame to rotate the sun's position vector
    R_MOD2J2000::SatelliteToolboxTransformations.DCM{TT} = r_eci_to_eci(
        MOD(), J2000(), time, eop_data
    )

    if body.name == :Sun
        pos_mod = sun_position_mod(time)
    elseif body.name == :Moon
        pos_mod = moon_position_mod(time)
    end

    return R_MOD2J2000 * pos_mod
end

"""
Computes the velocity of the celestial body using Vallado's ephemeris

# Arguments
- `ephem_type::Symbol`: Ephemeris type used to compute body's position.
- `body::CelestialBody`: Celestial body acting on the craft.
- `time::Number`: Current time of the simulation in seconds.

# Returns
- `body_position: SVector{3}`: The 3-dimensional third body position in the J2000 frame.
"""
function get_velocity(
    ephem_type::Vallado, body::CelestialBody, eop_data::T, time::TT
) where {T<:Union{EopIau1980,EopIau2000A,Nothing},TT}
    if body.name == :Sun
        vel_mod = sun_velocity_mod(time)
    else
        throw(
            ArgumentError("Vallado velocity ephemeris is only supported by Sun Currently")
        )
    end

    # Compute the MOD frame in the J2000 frame to rotate the sun's position vector
    R_MOD2J2000::SatelliteToolboxTransformations.DCM{TT} = r_eci_to_eci(
        MOD(), J2000(), time, eop_data
    )

    return R_MOD2J2000 * vel_mod
end

#TODO: ADD FULL STATE WITH SPICE SUPPORT
"""
Convenience to compute the ephemeris position of a CelestialBody in a ThirdBodyModel
Wraps get_position().

# Arguments
- `time::Number`: Current time of the simulation in seconds.

# Returns
- `body_position: SVector{3}`: The 3-dimensional third body position in the J2000 frame.

"""
function (model::ThirdBodyModel)(t::Number, return_type::Position)
    return get_position(model.ephem_type, model.body, model.eop_data, t)
end

function (model::ThirdBodyModel)(t::Number, return_type::Velocity)
    return get_velocity(model.ephem_type, model.body, model.eop_data, t)
end
