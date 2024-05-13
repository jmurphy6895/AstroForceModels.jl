# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Third Body Model and Ephemeris Functions
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export Vallado
Vallado() = :Vallado

export ThirdBodyModel
"""
Third Body Model Astro Model struct
Contains information to compute the acceleration of a third body force acting on a spacecraft.

# Fields
- `body::CelestialBody`: Celestial body acting on the craft.
- `ephem_type::Symbol`: Ephemeris type used to compute body's position. Options are currently :Vallado.
"""
@with_kw struct ThirdBodyModel{T} <: AbstractNonPotentialBasedForce where {
    T<:Union{EopIau1980,EopIau2000A,Nothing}
}
    body::CelestialBody
    eop_data::T = nothing
    ephem_type::Symbol = Vallado()
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
    ephem_type::Val{:Vallado}, body::CelestialBody, eop_data::T, time::Number
) where {T<:Union{EopIau1980,EopIau2000A,Nothing}}
    if body.name == :Sun
        pos_mod = sun_position_mod(time)
    elseif body.name == :Moon
        pos_mod = moon_position_mod(time)
    else
        throw(
            ArgumentError("Vallado ephemeris is only supported by Sun and Moon Currently")
        )
    end

    # Compute the MOD frame in the J2000 frame to rotate the sun's position vector
    R_MOD2J2000 = r_eci_to_eci(MOD(), J2000(), time, eop_data)

    return R_MOD2J2000 * pos_mod
end

@valsplit 1 function get_position(
    ephem_type::Symbol, body::CelestialBody, eop_data::T, time::Number
) where {T<:Union{EopIau1980,EopIau2000A,Nothing}}
    throw(ArgumentError("$ephem_type is not supported. Current options are :Vallado"))
end

"""
Convience to compute the ephemeris postion of a CelestialBody in a ThirdBodyModel
Wraps get_position().

# Arguments
- `time::Number`: Current time of the simulation in seconds.

# Returns
- `body_position: SVector{3}`: The 3-dimensional third body position in the J2000 frame.

"""
function (model::ThirdBodyModel)(t::Number)
    return get_position(model.ephem_type, model.body, model.eop_data, t)
end
