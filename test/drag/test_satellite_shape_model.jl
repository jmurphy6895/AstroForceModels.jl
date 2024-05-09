@testset "Cannonball Fixed Drag Satellite Shape" verbose = true begin
    radius = 10.0 * abs(randn())
    mass = 260.0 * abs(randn())
    drag_coeff = 2.0 * abs(randn())

    area = π * (radius^2.0)
    BC = drag_coeff * area / mass

    drag_model1 = CannonballFixedDrag(BC)
    drag_model2 = CannonballFixedDrag(radius, mass, drag_coeff)

    u = randn(6)
    p = [0.0]
    t = 1000.0 * abs(randn())

    @test ballistic_coefficient(u, p, t, drag_model1) == BC
    @test ballistic_coefficient(u, p, t, drag_model2) == BC

    @test drag_model1.radius == 1.0
    @test drag_model1.mass == 1.0
    @test drag_model1.drag_coeff == BC

    @test drag_model2.radius == radius
    @test drag_model2.mass == mass
    @test drag_model2.drag_coeff == drag_coeff
end

@testset "Cannonball Fixed Drag Satellite Shape Errors" verbose = true begin
    @test_throws ArgumentError CannonballFixedDrag(-1.0)
    @test_throws ArgumentError CannonballFixedDrag(-10.0, 260.0, 1.0)
    @test_throws ArgumentError CannonballFixedDrag(10.0, -260.0, 1.0)
    @test_throws ArgumentError CannonballFixedDrag(10.0, 260.0, -1.0)
end
