export R_SUN,
    R_EARTH, SPEED_OF_LIGHT, SOLAR_FLUX, μ_MOON, μ_SUN, EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS

# Radius of the Sun [km]
const R_SUN::Float64 = 6.955E5
# Radius of the Earth [km]
const R_EARTH::Float64 = 6378.1363
# Speed of Light [m/s]
const SPEED_OF_LIGHT::Float64 = 2.99792458E8
# Solar Flux
const SOLAR_FLUX::Float64 = 1360.8 / SPEED_OF_LIGHT
# Graviational Parameter of the Moon
const μ_MOON::Float64 = 4.902800118457551E3
# Graviational Parameter of the Sun
const μ_SUN::Float64 = 1.3271244004127946E11
# Earth's Angular Momentum per Unit Mass [kg*km^2/s]
const EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS::Float64 = 0.4 * R_EARTH^2 * EARTH_ANGULAR_SPEED
