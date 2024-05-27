# https://link.springer.com/article/10.1007/s10569-021-10014-y
export RelativityModel,
    relativity_accel,
    lense_thirring_accleration,
    schwartzchild_acceleration,
    de_Sitter_acceleration
#TODO: CAN THIS BE INCORPORATED INTO THE POTENTIAL FORCES IN DROMO
"""
Relativity Astro Model struct
Contains information to compute the acceleration of relativity acting on a spacecraft.

# Fields
- `central_body::ThirdBodyModel`: The data to compute the central body's graviational parameter.
- `sun_body::ThirdBodyModel`: The data to compute the Sun's position, velocity, and graviational parameter.
- `c::Number`: The speed of light [km/s].
- `γ::Number`: 
- `β::Number`:
- `schwartzchild::Bool`: Include the Schwartzchild relativity effect.
- `lense_thirring::Bool`: Include the Lense Thirring relativity effect.
- `de_Sitter::Bool`: Include the de Sitter relativity effect.
"""
@with_kw struct RelativityModel <: AbstractNonPotentialBasedForce
    central_body::ThirdBodyModel = ThirdBodyModel(;
        body=EarthBody(), eop_data=fetch_iers_eop()
    )
    sun_body::ThirdBodyModel = ThirdBodyModel(; body=SunBody(), eop_data=fetch_iers_eop())
    eop_data::Union{EopIau1980,EopIau2000A} = fetch_iers_eop()
    c::Number = SPEED_OF_LIGHT / 1E3
    γ::Number = 1.0
    β::Number = 1.0
    schwartzchild_effect::Bool = true
    lense_thirring_effect::Bool = true
    de_Sitter_effect::Bool = true
end

"""
    acceleration(u::AbstractArray, p::ComponentVector, t::Number, srp_model::SRPAstroModel)

Computes the srp acceleration acting on a spacecraft given a srp model and current state and 
parameters of an object.

# Arguments
- `u::AbstractArray`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `srp_model::SRPAstroModel`: SRP model struct containing the relevant information to compute the acceleration.

# Returns
- `acceleration: SVector{3}`: The 3-dimensional srp acceleration acting on the spacecraft.

"""
function acceleration(
    u::AbstractArray, p::ComponentVector, t::Number, relativity_model::RelativityModel
)
    J = SVector{3}(
        normalize(
            r_ecef_to_eci(ITRF(), J2000(), p.JD + t / 86400.0, relativity_model.eop_data)[
                :, 3
            ],
        ) * EARTH_ANGULAR_MOMENTUM_PER_UNIT_MASS,
    )

    sun_pos = relativity_model.sun_body(p.JD + t / 86400.0) ./ 1E3
    sun_vel = relativity_model.sun_body(p.JD + t / 86400.0; return_type=:velocity) ./ 1E3

    return relativity_accel(
        u,
        sun_pos,
        sun_vel,
        relativity_model.central_body.body.μ,
        relativity_model.sun_body.body.μ,
        J;
        c=relativity_model.c,
        γ=relativity_model.γ,
        β=relativity_model.β,
        schwartzchild_effect=relativity_model.schwartzchild_effect,
        lense_thirring_effect=relativity_model.lense_thirring_effect,
        de_Sitter_effect=relativity_model.de_Sitter_effect,
    )
end

"""


"""
function relativity_accel(
    u::AbstractArray,
    r_sun::AbstractArray,
    v_sun::AbstractArray,
    μ_body::Number,
    μ_Sun::Number,
    J::AbstractArray;
    c::Number=SPEED_OF_LIGHT / 1E3,
    γ::Number=1.0,
    β::Number=1.0,
    schwartzchild_effect::Bool=true,
    lense_thirring_effect::Bool=true,
    de_Sitter_effect::Bool=true,
)
    return SVector{3}(
        schwartzchild_effect * schwartzchild_acceleration(u, μ_body; c=c, γ=γ, β=β) +
        lense_thirring_effect * lense_thirring_acceleration(u, μ_body, J; c=c, γ=γ) +
        de_Sitter_effect * de_Sitter_acceleration(u, r_sun, v_sun, μ_Sun; c=c, γ=γ),
    )
end

"""


"""
function schwartzchild_acceleration(
    u::AbstractArray, μ_body::Number; c::Number=SPEED_OF_LIGHT, γ::Number=1.0, β::Number=1.0
)
    r = @view(u[1:3])
    r_norm = norm(r)
    ṙ = @view(u[4:6])

    schwartzchild = SVector{3}(
        μ_body / ((c^2.0) * (r_norm^3.0)) * (
            ((2.0 * (β + γ)) * (μ_body / r_norm) - γ * dot(ṙ, ṙ)) * r +
            2.0 * (1.0 + γ) * dot(r, ṙ) * ṙ
        ),
    )

    return schwartzchild
end

"""


"""
function lense_thirring_acceleration(
    u::AbstractArray,
    μ_body::Number,
    J::AbstractArray;
    c::Number=SPEED_OF_LIGHT,
    γ::Number=1.0,
)
    r = @view(u[1:3])
    r_norm = norm(r)
    ṙ = @view(u[4:6])

    lense_thirring = SVector{3}(
        (1.0 + γ) *
        (μ_body / ((c^2.0) * (r_norm^3.0))) *
        ((3.0 / r_norm^2) * cross(r, ṙ) * dot(r, J) + cross(ṙ, J)),
    )

    return lense_thirring
end

"""


"""
function de_Sitter_acceleration(
    u::AbstractArray,
    r_sun::AbstractArray,
    v_sun::AbstractArray,
    μ_Sun::Number;
    c::Number=SPEED_OF_LIGHT,
    γ::Number=1.0,
)
    ṙ = @view(u[4:6])

    de_sitter = SVector{3}(
        (1.0 + 2.0 * γ) *
        (-μ_Sun / ((c^2.0) * (norm(-r_sun)^3.0))) *
        cross(cross(-v_sun, -r_sun), ṙ),
    )

    return de_sitter
end
