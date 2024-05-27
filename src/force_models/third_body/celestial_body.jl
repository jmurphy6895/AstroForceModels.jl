#TODO: THIS PROBABLY BELONGS IN SATELLITE TOOLKIT CELESTIAL BODY
struct CelestialBody{T<:Number}
    name::Symbol
    central_body::Symbol
    jpl_code::Int
    μ::T
    Req::T
end

#TODO: NEED TO DETERMINE HOW TO GENERALIZE & ADD MORE BODIES
export CelestialBody
"""
Creates a CelestialBody Object from a name
Available: {Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, }

Arguments:
- `name::String`: Name of the Planet Object (Casing Doesn't Matter)
- (Optional) `T::DataType`: Datatype of PlanetBody Fields

Returns:
-`CelestialBody`: PlanetBody Object with Fields:
    - `name::String`: Name of Object
    - `central_body::String`: Name of Central Body Being Orbited
    - `jpl_code::Int`: NAIF ID Code
    - `μ::T`: Graviational Parameter [km/s]
    - `Req::T`: Equatorial Radius [km]
"""
function CelestialBody(name::Val{:Sun}; T::DataType=Float64)
    return CelestialBody{T}(
        :Sun,                           # Name
        :None,                          # Central Body
        1,                              # NAIF ID Code
        T(1.32712440018E11),            # μ [km^2/s^3]
        T(6.955E5),                     # Equatorial Radius [km]
    )
end

function CelestialBody(name::Val{:Earth}; T::DataType=Float64)
    return CelestialBody{T}(
        :Earth,                         # Name
        :Sun,                           # Central Body
        399,                            # NAIF ID Code
        T(3.986004415e5),               # μ [km^2/s^3]
        T(6378.1363),                   # Equatorial Radius [km]
    )
end

function CelestialBody(name::Val{:Moon}; T::DataType=Float64)
    return CelestialBody{T}(
        :Moon,                          # Name
        :Earth,                         # Central Body
        301,                            # NAIF ID Code
        T(4.90486959E3),                # μ [km^2/s^3]
        T(1738.1),                      # Equatorial Radius [km]
    )
end

@valsplit 1 function CelestialBody(name::Symbol; T::DataType=Float64)
    throw(
        ArgumentError(
            "Celestial body $name not recognized or not yet supported. Consider creating custom CelestialBody instead.",
        ),
    )
end

export SunBody, MoonBody
SunBody(; T::DataType=Float64) = CelestialBody(:Sun; T=T)
EarthBody(; T::DataType=Float64) = CelestialBody(:Earth; T=T)
MoonBody(; T::DataType=Float64) = CelestialBody(:Moon; T=T)
