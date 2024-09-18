#TODO: THIS PROBABLY BELONGS IN SATELLITE TOOLKIT CELESTIAL BODY
export CelestialBody
"""
Creates a CelestialBody Object from a name
Available: {Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, 

Fields:
-`CelestialBody{T<:Number}`: PlanetBody Object with Fields:
    - `name::String`: Name of Object
    - `central_body::String`: Name of Central Body Being Orbited
    - `jpl_code::Int`: NAIF ID Code
    - `μ::T`: Graviational Parameter [km/s]
    - `Req::T`: Equatorial Radius [km]
"""
struct CelestialBody{T<:Number}
    name::Symbol
    central_body::Symbol
    jpl_code::Int
    μ::T
    Req::T
end

export SunBody, MoonBody, EarthBody
function SunBody(; T::DataType=Float64)
    return CelestialBody{T}(
        :Sun,                           # Name
        :None,                          # Central Body
        1,                              # NAIF ID Code
        T(μ_SUN),                       # μ [km^2/s^3]
        T(R_SUN),                       # Equatorial Radius [km]
    )
end

function EarthBody(; T::DataType=Float64)
    return CelestialBody{T}(
        :Earth,                         # Name
        :Sun,                           # Central Body
        399,                            # NAIF ID Code
        T(μ_EARTH),                     # μ [km^2/s^3]
        T(R_EARTH),                     # Equatorial Radius [km]
    )
end

function MoonBody(; T::DataType=Float64)
    return CelestialBody{T}(
        :Moon,                          # Name
        :Earth,                         # Central Body
        301,                            # NAIF ID Code
        T(μ_MOON),                      # μ [km^2/s^3]
        T(R_MOON),                      # Equatorial Radius [km]
    )
end
