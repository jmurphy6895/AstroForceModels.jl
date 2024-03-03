#TODO: THIS PROBABLY BELONGS IN SATELLITE TOOLKIT CELESTIAL BODY
struct CelestialBody{T<:Number,V<:Number}
    name::String
    central_body::String
    jpl_code::Int64
    ephemeris::Function
    epoch::String

    μ::T
    Req::V
end

export CelestialBody
"""
Creates a CelestialBody Object from a name
Available: {Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, }

Arguments:
-'name::String': Name of the Planet Object (Casing Doesn't Matter)
-(Optional) 'T::DataType': Datatype of PlanetBody Fields

Returns:
-'CelestialBody': PlanetBody Object with Fields:
    name::String: Name of Object
    central_body::String: Name of Central Body Being Orbited
    jpl_code::Int64: NAIF ID Code
    ephemeris::Function: Ephemeris Object 

    μ::T: Graviational Parameter [km/s]
    Req::T: Equatorial Radius [km]

"""
function CelestialBody(name::String; T::DataType=Float64)
    if lowercase(name) == "sun"
        println("Creating Sun")
        return SunBody{T,T}(
            "Sun",                          # Name
            "",                             # Central Body
            1,                              # NAIF ID Code
            "",                             # Ephemeris
            1.32712440018E11,               # μ (km/s)
            1.0,
        )                            # Equatorial Radius (km)

    elseif lowercase(name) == "mercury"
        println("Creating Mercury")
        return PlanetBody{T,T}(
            "Mercury",                      # Name
            "Sun",                          # Central Body
            199,                            # NAIF ID Code
            "",                             # Ephemeris
            2.2032E4,                       # μ (km/s)
            2439.7,
        )                         # Equatorial Radius (km)

    elseif lowercase(name) == "venus"
        println("Creating Venus")
        return PlanetBody{T,T}(
            "Venus",                        # Name
            "Sun",                          # Central Body
            299,                            # NAIF ID Code
            "",                             # Ephemeris
            3.24858599E5,                   # μ (km/s)
            6051.8,
        )                         # Equatorial Radius (km)

    elseif lowercase(name) == "earth"
        println("Creating Earth")
        return PlanetBody{T,T}(
            "Earth",                        # Name
            "Sun",                          # Central Body
            399,                            # NAIF ID Code
            "",                             # Ephemeris
            3.98600433E5,                   # μ (km/s)
            6051.8,
        )                         # Equatorial Radius (km)

    elseif lowercase(name) == "mars"
        println("Creating Mars")
        return PlanetBody{T,T}(
            "Mars",                         # Name
            "Sun",                          # Central Body
            499,                            # NAIF ID Code
            "",                             # Ephemeris
            4.28283100E4,                   # μ (km/s)
            3396.19,
        )                        # Equatorial Radius (km)

    elseif lowercase(name) == "jupiter"
        println("Creating Jupiter")
        return PlanetBody{T,T}(
            "Jupiter",                      # Name
            "Sun",                          # Central Body
            599,                            # NAIF ID Code
            "",                             # Ephemeris
            1.26686536E8,                   # μ (km/s)
            71492.0,
        )                        # Equatorial Radius (km)

    elseif lowercase(name) == "saturn"
        println("Creating Saturn")
        return PlanetBody{T,T}(
            "Saturn",                       # Name
            "Sun",                          # Central Body
            699,                            # NAIF ID Code
            "",                             # Spice File
            3.7931208E7,                    # μ (km/s)
            60268.0,
        )                        # Equatorial Radius (km)

    elseif lowercase(name) == "uranus"
        println("Creating Uranus")
        return PlanetBody{T,T}(
            "Uranus",                       # Name
            "Sun",                          # Central Body
            799,                            # NAIF ID Code
            "",                             # Spice File
            5.7939513E6,                    # μ (km/s)
            25559.0,
        )                        # Equatorial Radius (km)

    elseif lowercase(name) == "neptune"
        println("Creating Neptune")
        return PlanetBody{T,T}(
            "Neptune",                      # Name
            "Sun",                          # Central Body
            899,                            # NAIF ID Code
            "",                             # Spice File
            6.835100E6,                     # μ (km/s)
            24764.0,
        )                        # Equatorial Radius (km)

    elseif lowercase(name) == "pluto"
        println("Creating Pluto")
        return PlanetBody{T,T}(
            "Pluto",                        # Name
            "Sun",                          # Central Body
            999,                            # NAIF ID Code
            "",                             # Spice File
            8.71E2,                         # μ (km/s)
            1188.3,
        )                         # Equatorial Radius (km)

    #TODO: UPDATE EPHEMERIS
    elseif lowercase(name) == "ceres"
        println("Creating Ceres")
        return PlanetBody{T,T}(
            "Ceres",                        # Name
            "Sun",                          # Central Body
            2000001,                        # NAIF ID Code
            "",                             # Spice File
            6.26325E1,                      # μ (km/s)
            476.2,
        )                          # Equatorial Radius (km)

    elseif lowercase(name) == "moon"
        println("Creating Moon")
        return CelestialBody{T,T}(
            "Moon",                         # Name
            "Earth",                        # Central Body
            301,                            # NAIF ID
            "de405.bsp",                    # Spice File
            "2001-01-01.5",                 # Epoch
            4.90486959E3,                   # μ
            0.0,                            # Req (km)
        )

    elseif lowercase(name) == "europa"
        println("Creating Europa")
        return CelestialBody{T,T}(
            "Europa",                       # Name
            "Jupiter",                      # Central Body
            502,                            # NAIF ID
            "jup365.bsp",                   # Spice File
            "2000-01-01.5",                 # Epoch
            3202.7,                         # μ
            0.0,                            # Req (km)
        )

    elseif lowercase(name) == "enceladus"
        println("Creating Enceladus")
        return CelestialBody{T,T}(
            "Enceladus",                    # Name
            "Saturn",                       # Central Body
            602,                            # NAIF ID
            "sat427.bsp",                   # Spice File
            "2000-01-01.5",                 # Epoch
            3202.7,                         # μ
            0.0,                            # Req (km)
        )

    else
        error(
            "Celestial body $name not recognized or not yet supported. Consider creating custom CelestialBody instead.",
        )
    end
end

export CustomCelestialBody
"""
Creates Custom CelestialBody Object from a name

Arguments:
-'name::String': Name of the Planet Object (Casing Doesn't Matter)
-'ephemeris::Function': Ephemeris function
-(Optional) 'central_body::String': Name of Central Body Being Orbited
-(Optional) 'jpl_code::Int64': NAIF ID Code
-(Optional) 'μ::Number': Graviational Parameter [km/s]
-(Optional) 'Req::Number': Equatorial Radius [km]

Returns:
-'CelestialBody': PlanetBody Object with Fields:
    name::String: Name of Object
    central_body::String: Name of Central Body Being Orbited
    jpl_code::Int64: NAIF ID Code
    ephemeris::Function: Ephemeris function
    μ::Number: Graviational Parameter [km/s]
    Req::Number: Equatorial Radius [km]
"""
function CustomCelestialBody(
    name::String,
    ephemeris::Function;
    central_body::String="",
    naif_id::Int=-1,
    μ::Number=0.0,
    Req::Number=0.0,
)
    print("Creating Custom Body: $name")
    return CelestialBody(name, central_body, naif_id, ephemeris, epoch, μ, Req)
end

function (body::CelestialBody)(t)
    return model.ephem(t)
end
