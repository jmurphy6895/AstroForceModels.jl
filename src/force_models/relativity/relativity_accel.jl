# https://link.springer.com/article/10.1007/s10569-021-10014-y

struct RelativityModel <: AbstractNonPotentialBasedForce
    central_body::ThirdBodyModel
    sun_body::ThirdBodyModel
end

function RelativityModel(central_body::Symbol)
    central_body = ThirdBodyModel(central_body)
    sun_body = ThirdBodyModel(:Sun)

    return RelativityModel(central_body, sun_body)
end

function acceleration(
    relativity_model::RelativityModel, u::AbstractArray, p::ComponentVector, t::Number
)

    #TODO: OFFER OPTION TO COMPUTE FROM EOP or SPICE EPHEMERIS 
    J = [0.0; 0.0; relativity_model.central_body.ω]

    sun_pos, sun_vel = relativity_model.sun_body.ephem(t)

    return relativity_accel(
        u,
        relativity_model.central_body.μ,
        relativity_model.sun_body.μ,
        sun_pos,
        sun_vel,
        p.c,
        J,
    )
end

"""


"""
function relativity_accel(
    u::AbstractArray,
    μ_body::Number,
    μ_Sun::Number,
    r_sun::AbstractArray,
    v_sun::AbstractArray,
    c::Number,
    J::AbstractArray;
    γ::Number=1.0,
    β::Number=1.0,
)
    r = @view(u[1:3])
    ṙ = @view(u[4:6])

    schwartzchild =
        μ_body / ((c^2.0) * (r_norm^3.0)) * (
            ((2.0 * (β + γ)) * (μ_body / r_norm) - γ * dot(ṙ, ṙ)) * r +
            2.0 * (1.0 + γ) * dot(r, ṙ) * ṙ
        )

    lense_thirring =
        (1.0 + γ) *
        (μ_body / ((c^2.0) * (r_norm^3.0))) *
        ((3.0 / r_norm) * cross(r, ṙ) * dot(r, J) + cross(ṙ, J))

    de_sitter =
        (1.0 + 2.0 * γ) *
        cross(cross(v_sun, (-μ_Sun / ((c^2.0) * (r_norm^3.0))) * r_sun), ṙ)

    return schwartzchild + lense_thirring + de_sitter
end
