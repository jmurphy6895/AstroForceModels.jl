export build_dynamics_model

"""
    build_dynamics_model(u::AbstractArray, p::ComponentVector, t::Number, models::AbstractArray{AbstractAstroForceModel})

Convenience funcction to compute the overall acceleration acting on the spacecraft.

# Arguments
- `u::AbstractArray`: Current State of the simulation.
- `p::ComponentVector`: Current parameters of the simulation.
- `t::Number`: Current time of the simulation.
- `models::AbstractArray{AbstractAstroForceModel}`: Array of acceleration models acting on the spacecraft

# Returns
- `acceleration: SVector{3}`: The 3-dimensional drag acceleration acting on the spacecraft.

"""
@generated function build_dynamics_model(u::AbstractArray, p::ComponentVector, t::Number, models::Vector{<:AbstractAstroForceModel})
    N = length(models.parameters)
    exprs = [:(acceleration(u, p, t, models[$i])) for i in 1:N]
    return :(SVector{3}($(foldl((a, b) -> :($a + $b), exprs))))
end