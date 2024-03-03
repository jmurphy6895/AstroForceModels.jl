struct ThirdBodyModel <: AbstractNonPotentialBasedForce
    body::CelestialBody
end

function ThirdBodyModel(name::Symbol)
    return ThirdBodyModel(CelestialBody(name))
end

function (model::ThirdBodyModel)(t::Number)
    return model.body(t)
end
