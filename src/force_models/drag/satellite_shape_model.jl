abstract type AbstractSatelliteDragModel end

struct Cannonball_Fixed_Drag <: AbstractSatelliteDragModel
    radius::Number
    mass::Number
    drag_coeff::Number
    ballistic_coeff::Number
end

function Cannonball_Fixed_Drag(ballistic_coeff::Number)
    return Cannonball_Fixed_Drag(1.0, 1.0, ballistic_coeff, ballistic_coeff)
end

function Cannonball_Fixed_Drag(radius::Number, mass::Number, drag_coeff::Number)

    area = π * radius^2.0

    return Cannonball_Fixed_Drag(radius, mass, drag_coeff, drag_coeff * area / mass)
end

function compute_reflectivity_coefficient(model::Cannonball_Fixed_Drag)
    return model.ballistic_coeff
end
