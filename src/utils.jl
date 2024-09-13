function angle_between_vectors(v1::AbstractVector, v2::AbstractVector)
    unitv1 = normalize(v1)
    unitv2 = normalize(v2)

    y = unitv1 .- unitv2
    x = unitv1 .+ unitv2

    a = 2 * atan(norm2(y), norm2(x))

    !(signbit(a) || signbit(float(T)(pi) - a)) ? a : (signbit(a) ? zero(T) : float(T)(pi))
end