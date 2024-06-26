# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Different Shadow Models used Mainly in SRP Calculation
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#TODO: REFERENCE
#   [1]
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
export Conical, Conical_Simplified, Cylindrical, INPE
Conical() = :Conical
Conical_Simplified() = :Conical_Simplified
Cylindrical() = :Cylindrical
INPE() = :INPE

export shadow_model
"""
    shadow_model(sat_pos::AbstractArray, sun_pos::AbstractArray, R_Sun::Number, R_Earth::Number, t::Number, ShadowModel::Symbol)
Computes the Lighting Factor of the Sun occur from the Umbra and Prenumbra of Earth's Shadow

# Arguments

- `sat_pos::AbstractArray`: The current satellite position.
- `sun_pos::AbstractArray`: The current Sun position.
- `R_Sun::Number`: The radius of the Sun.
- `R_Earth::Number`: The radius of the Earth.
- `ShadowModel::Symbol`: The Earth shadow model to use. Current Options -- :Cylindrical, :Conical, :Conical_Simplified, :None

# Returns

- `SVector{3}{Number}`: Inertial acceleration from the 3rd body
"""
@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    ShadowModel::Val{:Cylindrical};
    R_Sun::Number=R_SUN,
    R_Earth::Number=R_EARTH,
)
    sat_pos = @view(sat_pos[1:3])

    # Compute dot product between sun and satellite positions
    dp_sun_sat = dot(normalize(sun_pos), sat_pos)

    return if (
        dp_sun_sat > 0.0 || norm(sat_pos - dp_sun_sat * normalize(sun_pos)) > R_Earth
    )
        1.0
    else
        0.0
    end
end

@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    ShadowModel::Val{:INPE};
    R_Sun::Number=R_SUN,
    R_Earth::Number=R_EARTH,
)
    sun_pos_mag = norm(sun_pos)
    sat_pos_mag = norm(sat_pos)

    # Umbra
    Xu = R_Earth * sun_pos_mag / (R_Sun - R_Earth)
    Au = asin(R_Earth / Xu)

    # Prenumbra
    Xp = R_Earth * sun_pos_mag / (R_Sun + R_Earth)
    Ap = asin(R_Earth / Xp)

    # Projection of the Satellite Position onto the Sun Direction
    r_sat_i = dot(sat_pos, normalize(sun_pos)) * sun_pos / sun_pos_mag

    # Distance of the Umbral Cone and the Spacecraft
    δi = norm(sat_pos - r_sat_i)

    # Umbra Cone Terminator at Spacecraft Location
    ep_i = (Xu - sat_pos_mag) * tan(Au)

    # Prenumbra Cone Terminator at Spacecraft Location
    Kp_i = (Xp + sat_pos_mag) * tan(Ap)

    # SRP Scaling Factor
    F = 1.0

    if dot(sat_pos / sat_pos_mag, sun_pos / sun_pos_mag) < 0.0
        if δi < ep_i
            F = 0.0
        elseif δi < Kp_i
            F = 1.0 - (δi - ep_i) / (Kp_i - ep_i)
        end
    end

    return F
end

@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    ShadowModel::Val{:Conical_Simplified};
    R_Sun::Number=R_SUN,
    R_Earth::Number=R_EARTH,
)
    R_spacecraft_Sun = sat_pos - sun_pos

    con_a = asin(R_Sun / norm(R_spacecraft_Sun))
    con_b = asin(R_Earth / norm(sat_pos))
    #TODO: REPLACE WITH AngleBetweenVectors FORK
    con_c = AngleBetweenVectors.angle(
        normalize(@view(R_spacecraft_Sun[1:3])), normalize(@view(sat_pos[1:3]))
    )

    return if con_c ≥ (con_b + con_a)
        1.0
    elseif con_c < (con_b - con_a)
        0.0
    else
        0.5 + (con_c - con_b) / (2.0 * con_a)
    end
end

@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    ShadowModel::Val{:Conical};
    R_Sun::Number=R_SUN,
    R_Earth::Number=R_EARTH,
)
    R_spacecraft_Sun = sat_pos - sun_pos

    con_a = asin(R_Sun / norm(R_spacecraft_Sun))
    con_b = asin(R_Earth / norm(sat_pos))
    #TODO: REPLACE WITH AngleBetweenVectors FORK
    con_c = AngleBetweenVectors.angle(
        normalize(@view(R_spacecraft_Sun[1:3])), normalize(@view(sat_pos[1:3]))
    )

    return if con_c ≥ (con_b + con_a)
        1.0
    elseif con_c < (con_b - con_a)
        0.0
    elseif con_c < (con_a - con_b)
        1.0 - (con_b^2) / (con_a^2)
    else
        x = (con_c^2 + con_a^2 - con_b^2) / (2.0 * con_c)
        y = √(con_a^2 - x^2)
        area = con_a^2 * acos(x / con_a) + con_b^2 * acos((con_c - x) / con_b) - con_c * y
        1.0 - area / (π * con_a^2)
    end
end

@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    ShadowModel::Val{:None};
    R_Sun::Number=R_SUN,
    R_Earth::Number=R_EARTH,
)
    return 1.0
end

@valsplit 3 function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    ShadowModel::Symbol;
    R_Sun::Number=R_SUN,
    R_Earth::Number=R_EARTH,
)
    return error(ArgumentError("Shadow Model is not Defined $ShadowModel"))
end
