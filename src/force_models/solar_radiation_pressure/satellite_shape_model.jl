abstract type AbstractSatelliteSRPModel end

struct Cannonball_Fixed_SRP <: AbstractSatelliteSRPModel
    radius::Number
    mass::Number
    reflectivity_coeff::Number
    reflectivity_ballistic_coeff::Number
end

function Cannonball_Fixed_SRP(reflectivity_ballistic_coeff::Number)
    return Cannonball_Fixed_SRP(
        1.0, 1.0, reflectivity_ballistic_coeff, reflectivity_ballistic_coeff
    )
end

function Cannonball_Fixed_SRP(radius::Number, mass::Number, reflectivity_coeff::Number)
    area = π * radius^2.0

    return Cannonball_Fixed_SRP(
        radius, mass, reflectivity_coeff, reflectivity_coeff * area / mass
    )
end

function compute_ballistic_coefficient(model::Cannonball_Fixed_SRP)
    return model.reflectivity_ballistic_coeff
end
